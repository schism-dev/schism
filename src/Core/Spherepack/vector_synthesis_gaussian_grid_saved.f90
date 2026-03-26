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
submodule(vector_synthesis_routines) vector_synthesis_gaussian_saved

contains

    ! Purpose:
    !
    !     subroutine vhsgs(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
    !                      mdab, ndab, wvhsgs, ierror)
    !
    !     subroutine vhsgs performs the vector spherical harmonic synthesis
    !     of the arrays br, bi, cr, and ci and stores the result in the
    !     arrays v and w.  the synthesis is performed on an equally spaced
    !     longitude grid and a gaussian colatitude grid (measured from
    !     the north pole). v(i, j) and w(i, j) are the colatitudinal and
    !     east longitudinal components respectively, located at the i(th)
    !     colatitude gaussian point (see nlat below) and longitude
    !     phi(j) = (j-1)*2*pi/nlon.  the spectral respresentation of (v, w)
    !     is given below at output parameters v, w.
    !
    !     input parameters
    !
    !     nlat   the number of points in the gaussian colatitude grid on the
    !            full sphere. these lie in the interval (0, pi) and are computed
    !            in radians in theta(1) <...< theta(nlat) by subroutine
    !            compute_gaussian_latitudes_and_weights.
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
    !     nt     the number of syntheses.  in the program that calls vhsgs, 
    !            the arrays v, w, br, bi, cr, and ci can be three dimensional
    !            in which case multiple syntheses will be performed.
    !            the third index is the synthesis index which assumes the
    !            values k=1, ..., nt.  for a single synthesis set nt=1. the
    !            discription of the remaining parameters is simplified
    !            by assuming that nt=1 or that all the arrays are two
    !            dimensional.
    !
    !     idvw   the first dimension of the arrays v, w as it appears in
    !            the program that calls vhags. if ityp <= 2 then idvw
    !            must be at least nlat.  if ityp > 2 and nlat is
    !            even then idvw must be at least nlat/2. if ityp > 2
    !            and nlat is odd then idvw must be at least (nlat + 1)/2.
    !
    !     jdvw   the second dimension of the arrays v, w as it appears in
    !            the program that calls vhsgs. jdvw must be at least nlon.
    !
    !     br, bi  two or three dimensional arrays (see input parameter nt)
    !     cr, ci  that contain the vector spherical harmonic coefficients
    !            in the spectral representation of v(i, j) and w(i, j) given
    !            below at the discription of output parameters v and w.
    !
    !     mdab   the first dimension of the arrays br, bi, cr, and ci as it
    !            appears in the program that calls vhsgs. mdab must be at
    !            least min(nlat, nlon/2) if nlon is even or at least
    !            min(nlat, (nlon + 1)/2) if nlon is odd.
    !
    !     ndab   the second dimension of the arrays br, bi, cr, and ci as it
    !            appears in the program that calls vhsgs. ndab must be at
    !            least nlat.
    !
    !     wvhsgs an array which must be initialized by subroutine vhsgsi.
    !            once initialized, wvhsgs can be used repeatedly by vhsgs
    !            as long as nlon and nlat remain unchanged.  wvhsgs must
    !            not be altered between calls of vhsgs.
    !
    !     lvhsgs the dimension of the array wvhsgs as it appears in the
    !            program that calls vhsgs. define
    !
    !               l1 = min(nlat, nlon/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lvhsgs must be at least
    !
    !                 l1*l2*(2*nlat-l1+1)+nlon+15+2*nlat
    !
    !     output parameters
    !
    !     v, w    two or three dimensional arrays (see input parameter nt)
    !            in which the synthesis is stored. v is the colatitudinal
    !            component and w is the east longitudinal component.
    !            v(i, j), w(i, j) contain the components at the guassian colatitude
    !            point theta(i) and longitude phi(j) = (j-1)*2*pi/nlon.
    !            the index ranges are defined above at the input parameter
    !            ityp. v and w are computed from the formulas given below.
    !
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
    !     7.   theta(i) = i(th) gaussian grid point and phi(j) = (j-1)*2*pi/nlon
    !
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
    !            = 9  error in the specification of lvhsgs
    !
    module subroutine vhsgs(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
        mdab, ndab, wvhsgs, ierror)

        ! Dummy arguments
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: nlon
        integer(ip), intent(in)     :: ityp
        integer(ip), intent(in)     :: nt
        real(wp),    intent(out)    :: v(idvw, jdvw, nt)
        real(wp),    intent(out)    :: w(idvw, jdvw, nt)
        integer(ip), intent(in)     :: idvw
        integer(ip), intent(in)     :: jdvw
        real(wp),    intent(in)     :: br(mdab, ndab, nt)
        real(wp),    intent(in)     :: bi(mdab, ndab, nt)
        real(wp),    intent(in)     :: cr(mdab, ndab, nt)
        real(wp),    intent(in)     :: ci(mdab, ndab, nt)
        integer(ip), intent(in)     :: mdab
        integer(ip), intent(in)     :: ndab
        real(wp),    intent(in)     :: wvhsgs(:)
        integer(ip), intent(out)    :: ierror

        ! Local variables
        integer(ip) :: idv, idz, imid, ist, lnl, lzimn, mmax
        integer(ip) :: lwork, workspace_indices(7)

        associate (lvhsgs => size(wvhsgs))

            imid = (nlat + 1)/2
            mmax = min(nlat, (nlon + 1)/2)
            idz = (mmax*(2*nlat-mmax+1))/2
            lzimn = idz*imid

            ! Check calling arguments
            if (nlat < 3) then
                ierror = 1
            else if (nlon < 1) then
                ierror = 2
            else if (ityp < 0 .or. ityp > 8) then
                ierror = 3
            else if (nt < 0) then
                ierror = 4
            else if ( &
                (ityp <= 2 .and. idvw < nlat) &
                .or. &
                (ityp > 2 .and. idvw < imid)) then
                ierror = 5
            else if (jdvw < nlon) then
                ierror = 6
            else if (mdab < mmax) then
                ierror = 7
            else if (ndab < nlat) then
                ierror = 8
            else if (lvhsgs < 2*lzimn+nlon+15) then
                ierror = 9
            else
                ierror = 0
            end if

            ! Check error flag
            if (ierror /= 0) return

            ! Set required workspace size
            select case (ityp)
                case(0:2)
                    ist = imid
                    idv = nlat
                case default
                    ist = 0
                    idv = imid
            end select

            lnl = nt*idv*nlon
            lwork = (2 * lnl) + (idv * nlon)

            block
                real(wp) :: work(lwork)

                !  Compute workspace indices
                workspace_indices = get_vhsgs_workspace_indices(nlat, imid, ist, lnl)

                associate (&
                    jw1 => workspace_indices(1), &
                    jw2 => workspace_indices(2), &
                    jw3 => workspace_indices(3), &
                    iw1 => workspace_indices(4), &
                    iw2 => workspace_indices(5), &
                    iw3 => workspace_indices(6), &
                    iw4 => workspace_indices(7) &
                    )
                    call vhsgs_lower_utility_routine(nlat, nlon, ityp, nt, imid, idvw, &
                        jdvw, v, w, mdab, ndab, br, bi, cr, ci, idv, work, work(iw1:), &
                        work(iw2:), work(iw3:), work(iw4:), idz, wvhsgs(jw1:), &
                        wvhsgs(jw2:), wvhsgs(jw3:))
                end associate
            end block
        end associate

    end subroutine vhsgs

    ! Purpose:
    !
    !     Computes the gaussian points theta, gauss
    !     weights wts, and the components vb and wb of the vector
    !     harmonics.
    !
    !     set imid = (nlat + 1)/2 and lmn=(nlat*(nlat + 1))/2 then
    !     wvhsgs must have 2*(imid*lmn+nlat)+nlon+15 locations
    !
    !     real array dwork must have
    !     3*nlat*(nlat + 1)+5*nlat+1 = nlat*(3*nlat+8)+1
    !     locations which is determined by the size of dthet, 
    !     dwts, dwork, and dpbar in vhsgs_lower_utility_routine
    !
    !     subroutine vhsgsi(nlat, nlon, wvhsgs, ierror)
    !
    !     subroutine vhsgsi initializes the array wvhsgs which can then be
    !     used repeatedly by subroutine vhsgs until nlat or nlon is changed.
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
    !     lvhsgs the dimension of the array wvhsgs as it appears in the
    !            program that calls vhsgs. define
    !
    !               l1 = min(nlat, nlon/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lvhsgs must be at least
    !
    !                 l1*l2*(2*nlat-l1+1)+nlon+15+2*nlat
    !
    !     output parameters
    !
    !     wvhsgs an array which is initialized for use by subroutine vhsgs.
    !            once initialized, wvhsgs can be used repeatedly by vhsgs
    !            as long as nlat and nlon remain unchanged.  wvhsgs must not
    !            be altered between calls of vhsgs.
    !
    !
    !     ierror = 0  no errors
    !            = 1  error in the specification of nlat
    !            = 2  error in the specification of nlon
    !            = 3  error in the specification of lvhsgs
    !
    module subroutine vhsgsi(nlat, nlon, wvhsgs, ierror)

        ! Dummy arguments
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: nlon
        real(wp),    intent(out)    :: wvhsgs(:)
        integer(ip), intent(out)    :: ierror

        ! Local variables
        integer(ip) :: imid, lmn, ldwork
        type(SpherepackUtility) :: util

        associate (lvhsgs => size(wvhsgs))

            imid = (nlat + 1)/2
            lmn = (nlat*(nlat + 1))/2

            !  Check calling arguments
            if (nlat < 3) then
                ierror = 1
            else if (nlon < 1) then
                ierror = 2
            else if (lvhsgs < 2*(imid*lmn) + nlon + 15) then
                ierror = 3
            else
                ierror = 0
            end if

            ! Check error flag
            if (ierror /= 0) return

            ! Set required workspace size
            ldwork = (nlat*3*(nlat+3)+2)/2

            block
                integer(ip) :: jw1, jw2, jw3
                integer(ip) :: iw1, iw2, iw3, iw4
                real(wp)    :: dwork(ldwork)

                ! Set saved workspace pointer indices
                jw1 = 1
                jw2 = jw1+imid*lmn
                jw3 = jw2+imid*lmn

                ! Set unsaved workspace pointer indices
                iw1 = 1
                iw2 = iw1+nlat
                iw3 = iw2+nlat
                iw4 = iw3+3*imid*nlat

                call vhgsi_lower_utility_routine(nlat, imid, wvhsgs(jw1:), &
                    wvhsgs(jw2:), dwork(iw1:), dwork(iw2:), dwork(iw3:), dwork(iw4:))

                call util%hfft%initialize(nlon, wvhsgs(jw3:))
            end block
        end associate

    end subroutine vhsgsi

    pure function get_vhsgs_workspace_indices(nlat, imid, ist, lnl) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: imid
        integer(ip), intent(in)  :: ist
        integer(ip), intent(in)  :: lnl
        integer(ip)              :: return_value(7)

        ! Local variables
        integer(ip) :: lmn

        associate (i => return_value)

            !  set wvhsgs pointers
            lmn = nlat*(nlat + 1)/2
            i(1) = 1
            i(2) = i(1)+imid*lmn
            i(3) = i(2)+imid*lmn

            !  set work pointers
            i(4) = ist+1
            i(5) = lnl+1
            i(6) = i(5)+ist
            i(7) = i(5)+lnl
        end associate

    end function get_vhsgs_workspace_indices

    subroutine vhsgs_lower_utility_routine(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
        ndab, br, bi, cr, ci, idv, ve, vo, we, wo, work, idz, vb, wb, wrfft)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: ityp
        integer(ip), intent(in)  :: nt
        integer(ip), intent(in)  :: imid
        integer(ip), intent(in)  :: idvw
        integer(ip), intent(in)  :: jdvw
        real(wp),    intent(out) :: v(idvw, jdvw, nt)
        real(wp),    intent(out) :: w(idvw, jdvw, nt)
        integer(ip), intent(in)  :: mdab
        integer(ip), intent(in)  :: ndab
        real(wp),    intent(in)  :: br(mdab, ndab, nt)
        real(wp),    intent(in)  :: bi(mdab, ndab, nt)
        real(wp),    intent(in)  :: cr(mdab, ndab, nt)
        real(wp),    intent(in)  :: ci(mdab, ndab, nt)
        integer(ip), intent(in)  :: idv
        real(wp),    intent(out) :: ve(idv, nlon, nt)
        real(wp),    intent(out) :: vo(idv, nlon, nt)
        real(wp),    intent(out) :: we(idv, nlon, nt)
        real(wp),    intent(out) :: wo(idv, nlon, nt)
        real(wp),    intent(out) :: work(*)
        integer(ip), intent(in)  :: idz
        real(wp),    intent(in)  :: vb(imid, *)
        real(wp),    intent(in)  :: wb(imid, *)
        real(wp),    intent(in)  :: wrfft(:)

        ! Local variables
        integer(ip)    :: i, imm1, k, m, mb, mmax
        integer(ip)    :: mn, mp1, mp2, odd_stride, even_stride, np1
        type(VectorSynthesisUtility) :: util

        call util%synthesis_setup(even_stride, imid, imm1, mmax, nlat, &
            odd_stride, ve, vo, we, wo)

        vector_symmetry_cases: select case (ityp)
            case (0)
                                      ! case ityp=0   no symmetries
                !
                ! case m = 0
                !
                do k=1, nt
                    do np1=2, even_stride, 2
                        do i=1, imid
                            ve(i, 1, k)=ve(i, 1, k)+br(1, np1, k)*vb(i, np1)
                            we(i, 1, k)=we(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do

                do k=1, nt
                    do np1=3, odd_stride, 2
                        do i=1, imm1
                            vo(i, 1, k)=vo(i, 1, k)+br(1, np1, k)*vb(i, np1)
                            wo(i, 1, k)=wo(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do

                ! case m = 1 through nlat-1
                if (2 <= mmax) then
                    do mp1=2, mmax
                        m = mp1-1
                        mb = m * (nlat-1) - (m * (m - 1))/2
                        mp2 = mp1+1
                        if (mp1 <= odd_stride) then
                            do k=1, nt
                                do np1=mp1, odd_stride, 2
                                    mn = mb+np1
                                    do i=1, imm1
                                        vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                        ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                        vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                        ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                        wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                        we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                        wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                        we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                    end do

                                    if (mod(nlat, 2) /= 0) then
                                        ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                            -ci(mp1, np1, k)*wb(imid, mn)
                                        ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                            +cr(mp1, np1, k)*wb(imid, mn)
                                        we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                            -bi(mp1, np1, k)*wb(imid, mn)
                                        we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                            +br(mp1, np1, k)*wb(imid, mn)
                                    end if
                                end do
                            end do
                        end if

                        if (mp2 <= even_stride) then
                            do k=1, nt
                                do np1=mp2, even_stride, 2
                                    mn = mb+np1
                                    do i=1, imm1
                                        ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                        vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                        ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                        vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                        we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                        wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                        we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                        wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                    end do
                                    if (mod(nlat, 2) /= 0) then
                                        ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                            +br(mp1, np1, k)*vb(imid, mn)
                                        ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                            +bi(mp1, np1, k)*vb(imid, mn)
                                        we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                            -cr(mp1, np1, k)*vb(imid, mn)
                                        we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                            -ci(mp1, np1, k)*vb(imid, mn)
                                    end if
                                end do
                            end do
                        end if
                    end do
                end if
            case (1)
                !
                ! case ityp=1   no symmetries,  cr and ci equal zero
                !
                ! case m = 0
                !
                do k=1, nt
                    do np1=2, even_stride, 2
                        do i=1, imid
                            ve(i, 1, k)=ve(i, 1, k)+br(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do

                do k=1, nt
                    do np1=3, odd_stride, 2
                        do i=1, imm1
                            vo(i, 1, k)=vo(i, 1, k)+br(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) exit vector_symmetry_cases
                do mp1=2, mmax
                    m = mp1-1
                    !     mb = m*(nlat-1)-(m*(m-1))/2
                    mb = m*nlat-(m*(m+1))/2
                    mp2 = mp1+1
                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do np1=mp1, odd_stride, 2
                                mn = mb+np1
                                do i=1, imm1
                                    vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                    vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                    we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                    we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                end do
                                if (mod(nlat, 2) /= 0) then
                                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                        -bi(mp1, np1, k)*wb(imid, mn)
                                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                        +br(mp1, np1, k)*wb(imid, mn)
                                end if
                            end do
                        end do
                    end if

                    if (mp2 <= even_stride) then
                        do k=1, nt
                            do np1=mp2, even_stride, 2
                                mn = mb+np1
                                do i=1, imm1
                                    ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                    ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                    wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                    wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                end do
                                if (mod(nlat, 2) /= 0) then
                                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                        +br(mp1, np1, k)*vb(imid, mn)
                                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                        +bi(mp1, np1, k)*vb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                end do
            case (2)
                !
                ! case ityp=2   no symmetries,  br and bi are equal to zero
                !
                ! case m = 0
                !
                do k=1, nt
                    do np1=2, even_stride, 2
                        do i=1, imid
                            we(i, 1, k)=we(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do

                do k=1, nt
                    do np1=3, odd_stride, 2
                        do i=1, imm1
                            wo(i, 1, k)=wo(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) exit vector_symmetry_cases

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*nlat-(m*(m+1))/2
                    mp2 = mp1+1

                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do np1=mp1, odd_stride, 2
                                mn = mb+np1
                                do i=1, imm1
                                    ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                    ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                    wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                    wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                end do
                                if (mod(nlat, 2) /= 0) then
                                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                        -ci(mp1, np1, k)*wb(imid, mn)
                                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                        +cr(mp1, np1, k)*wb(imid, mn)
                                end if
                            end do
                        end do
                    end if

                    if (mp2 <= even_stride) then
                        do k=1, nt
                            do np1=mp2, even_stride, 2
                                mn = mb+np1
                                do i=1, imm1
                                    vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                    vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                    we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                    we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                end do
                                if (mod(nlat, 2) /= 0) then
                                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                        -cr(mp1, np1, k)*vb(imid, mn)
                                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                        -ci(mp1, np1, k)*vb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                end do
            case (3)
                !
                ! case ityp=3   v even,  w odd
                !
                ! case m = 0
                !
                do k=1, nt
                    do np1=2, even_stride, 2
                        do i=1, imid
                            ve(i, 1, k)=ve(i, 1, k)+br(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do

                do k=1, nt
                    do np1=3, odd_stride, 2
                        do i=1, imm1
                            wo(i, 1, k)=wo(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) exit vector_symmetry_cases

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*nlat-(m*(m+1))/2
                    mp2 = mp1+1
                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do np1=mp1, odd_stride, 2
                                mn = mb+np1
                                do i=1, imm1
                                    ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                    ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                    wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                    wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                end do
                                if (mod(nlat, 2) /= 0) then
                                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                        -ci(mp1, np1, k)*wb(imid, mn)
                                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                        +cr(mp1, np1, k)*wb(imid, mn)
                                end if
                            end do
                        end do
                    end if

                    if (mp2 <= even_stride) then
                        do k=1, nt
                            do np1=mp2, even_stride, 2
                                mn = mb+np1
                                do i=1, imm1
                                    ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                    ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                    wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                    wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                end do
                                if (mod(nlat, 2) /= 0) then
                                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                        +br(mp1, np1, k)*vb(imid, mn)
                                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                        +bi(mp1, np1, k)*vb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                end do
            case (4)
                !
                ! case ityp=4   v even,  w odd, and both cr and ci equal zero
                !
                ! case m = 0
                !
                do k=1, nt
                    do np1=2, even_stride, 2
                        do i=1, imid
                            ve(i, 1, k)=ve(i, 1, k)+br(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) exit vector_symmetry_cases

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*nlat-(m*(m+1))/2
                    mp2 = mp1+1

                    if (mp2 <= even_stride) then
                        do k=1, nt
                            do np1=mp2, even_stride, 2
                                mn = mb+np1
                                do i=1, imm1
                                    ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                    ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                    wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                    wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                end do
                                if (mod(nlat, 2) /= 0) then
                                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                        +br(mp1, np1, k)*vb(imid, mn)
                                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                        +bi(mp1, np1, k)*vb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                end do
            case (5)
                 !
                 ! case ityp=5   v even,  w odd,     br and bi equal zero
                 !
                 ! case m = 0
                 !
                do k=1, nt
                    do np1=3, odd_stride, 2
                        do i=1, imm1
                            wo(i, 1, k)=wo(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) exit vector_symmetry_cases

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*nlat-(m*(m+1))/2
                    mp2 = mp1+1

                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do np1=mp1, odd_stride, 2
                                mn = mb+np1
                                do i=1, imm1
                                    ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                    ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                    wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                    wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                end do
                                if (mod(nlat, 2) /= 0) then
                                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                        -ci(mp1, np1, k)*wb(imid, mn)
                                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                        +cr(mp1, np1, k)*wb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                end do
            case (6)
                !
                ! case ityp=6   v odd  ,  w even
                !
                ! case m = 0
                !
                do k=1, nt
                    do np1=2, even_stride, 2
                        do i=1, imid
                            we(i, 1, k)=we(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do

                do k=1, nt
                    do np1=3, odd_stride, 2
                        do i=1, imm1
                            vo(i, 1, k)=vo(i, 1, k)+br(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) exit vector_symmetry_cases

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*nlat-(m*(m+1))/2
                    mp2 = mp1+1
                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do np1=mp1, odd_stride, 2
                                mn = mb+np1
                                do i=1, imm1
                                    vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                    vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                    we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                    we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                end do
                                if (mod(nlat, 2) /= 0) then
                                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                        -bi(mp1, np1, k)*wb(imid, mn)
                                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                        +br(mp1, np1, k)*wb(imid, mn)
                                end if
                            end do
                        end do
                    end if

                    if (mp2 <= even_stride) then
                        do k=1, nt
                            do np1=mp2, even_stride, 2
                                mn = mb+np1
                                do i=1, imm1
                                    vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                    vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                    we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                    we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                end do
                                if (mod(nlat, 2) /= 0) then
                                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                        -cr(mp1, np1, k)*vb(imid, mn)
                                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                        -ci(mp1, np1, k)*vb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                end do
            case (7)
                !
                ! case ityp=7   v odd, w even   cr and ci equal zero
                !
                ! case m = 0
                !
                do k=1, nt
                    do np1=3, odd_stride, 2
                        do i=1, imm1
                            vo(i, 1, k)=vo(i, 1, k)+br(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) exit vector_symmetry_cases

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*nlat-(m*(m+1))/2
                    mp2 = mp1+1
                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do np1=mp1, odd_stride, 2
                                mn = mb+np1
                                do i=1, imm1
                                    vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                    vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                    we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                    we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                end do
                                if (mod(nlat, 2) /= 0) then
                                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                        -bi(mp1, np1, k)*wb(imid, mn)
                                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                        +br(mp1, np1, k)*wb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                end do
            case (8)
                !
                ! case ityp=8   v odd,  w even   br and bi equal zero
                !
                ! case m = 0
                !
                do k=1, nt
                    do np1=2, even_stride, 2
                        do i=1, imid
                            we(i, 1, k)=we(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) exit vector_symmetry_cases

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*nlat-(m*(m+1))/2
                    mp2 = mp1+1
                    if (mp2 <= even_stride) then
                        do k=1, nt
                            do np1=mp2, even_stride, 2
                                mn = mb+np1
                                do i=1, imm1
                                    vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                    vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                    we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                    we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                end do
                                if (mod(nlat, 2) /= 0) then
                                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                        -cr(mp1, np1, k)*vb(imid, mn)
                                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                        -ci(mp1, np1, k)*vb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                end do
        end select vector_symmetry_cases

        call util%assemble_transform(idvw, jdvw, idv, imid, &
            imm1, ityp, nlat, nlon, nt, v, ve, vo, w, we, wo, wrfft)

    end subroutine vhsgs_lower_utility_routine

    subroutine vhgsi_lower_utility_routine(nlat, imid, vb, wb, dthet, dwts, dpbar, work)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: imid
        real(wp),    intent(out) :: vb(imid, *)
        real(wp),    intent(out) :: wb(imid, *)
        real(wp),    intent(out) :: dthet(nlat)
        real(wp),    intent(out) :: dwts(nlat)
        real(wp),    intent(out) :: dpbar(imid, nlat, 3)
        real(wp),    intent(out) :: work(:)

        ! Local variables
        integer(ip)         :: i, id, ierror, ix, iy
        integer(ip)         :: m, mn, n, nm, np, nz
        real(wp)            :: abel, bbel, cbel, dcf
        type(SpherepackUtility) :: util

        !  Compute gauss points and weights
        call compute_gaussian_latitudes_and_weights(nlat, dthet, dwts, ierror)
        !
        !  Compute associated legendre functions
        !
        !    Compute m=n=0 legendre polynomials for all theta(i)
        !
        dpbar(:, 1, 1) = cos(PI/4)
        vb(:, 1) = ZERO
        wb(:, 1) = ZERO
        !
        !  main loop for remaining vb, and wb
        !
        do n=1, nlat-1
            nm = mod(n - 2, 3) + 1
            nz = mod(n - 1, 3) + 1
            np = mod(n, 3) + 1
            !
            !  Compute dpbar for m=0
            !
            call util%compute_fourier_coefficients(0, n, work)
            mn = get_index(0, n, nlat)
            do i=1, imid
                call util%compute_legendre_polys_from_fourier_coeff(0, n, dthet(i), work, dpbar(i, 1, np))
            end do
            !
            !  Compute dpbar for m=1
            !
            call util%compute_fourier_coefficients(1, n, work)
            mn = get_index(1, n, nlat)
            do i=1, imid
                call util%compute_legendre_polys_from_fourier_coeff(1, n, dthet(i), work, dpbar(i, 2, np))
            end do
            !
            !  Compute and store dpbar for m=2, n
            !
            if (2 <= n) then
                do m=2, n
                    abel = sqrt(real((2*n+1)*(m+n-2)*(m+n-3), kind=wp)/ &
                        real((2*n-3)*(m+n-1)*(m+n), kind=wp))
                    bbel = sqrt(real((2*n+1)*(n-m-1)*(n-m), kind=wp)/ &
                        real((2*n-3)*(m+n-1)*(m+n), kind=wp))
                    cbel = sqrt(real((n-m+1)*(n-m+2), kind=wp)/ &
                        real((m+n-1)*(m+n), kind=wp))
                    id = get_index(m, n, nlat)
                    if (m < n-1) then
                        dpbar(:, m+1, np) = abel*dpbar(:, m-1, nm)+bbel*dpbar(:, m+1, nm)-cbel*dpbar(:, m-1, np)
                    else
                        dpbar(:, m+1, np) = abel*dpbar(:, m-1, nm)-cbel*dpbar(:, m-1, np)
                    end if
                end do
            end if
            !
            !  Compute the derivative of the functions
            !
            ix = get_index(0, n, nlat)
            iy = get_index(n, n, nlat)
            do i=1, imid
                vb(i, ix) = -dpbar(i, 2, np)
                vb(i, iy) = dpbar(i, n, np)/sqrt(real(2*(n+1), kind=wp))
            end do
            !
            if (n /= 1) then
                dcf = sqrt(real(4*n*(n+1), kind=wp))
                do m=1, n-1
                    ix = get_index(m, n, nlat)
                    abel = sqrt(real((n+m)*(n-m+1), kind=wp))/dcf
                    bbel = sqrt(real((n-m)*(n+m+1), kind=wp))/dcf
                    do i=1, imid
                        vb(i, ix) = abel*dpbar(i, m, np)-bbel*dpbar(i, m+2, np)
                    end do
                end do
            end if
            !
            !     compute the vector harmonic w(theta) = m*pbar/cos(theta)
            !
            !     set wb=0 for m=0
            !
            ix = get_index(0, n, nlat)
            wb(1:imid, ix) = ZERO
            !
            !     compute wb for m=1, n
            !
            dcf = sqrt(real(2*n+1, kind=wp)/real(4*n*(n+1)*(n+n-1), kind=wp))

            do m=1, n
                ix = get_index(m, n, nlat)
                abel = dcf*sqrt(real((n+m)*(n+m-1), kind=wp))
                bbel = dcf*sqrt(real((n-m)*(n-m-1), kind=wp))
                if (m < n-1) then
                    wb(:, ix) = abel*dpbar(:, m, nz) + bbel*dpbar(:, m+2, nz)
                else
                    wb(:, ix) = abel*dpbar(:, m, nz)
                end if
            end do
        end do

    end subroutine vhgsi_lower_utility_routine

    pure function get_index(m, n, nlat) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: m
        integer(ip), intent(in) :: n
        integer(ip), intent(in) :: nlat
        integer(ip)             :: return_value

        return_value = m*nlat-(m*(m+1))/2+n+1

    end function get_index

end submodule vector_synthesis_gaussian_saved
