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
!
! ... file sfvpgc.f
!
!     this file includes documentation and code for
!     subroutine sfvpgc          i
!
! ... files which must be loaded with sfvpgc.f
!
!     type_SpherepackUtility.f, type_RealPeriodicFastFourierTransform.f, vhagc.f, shsgc.f, compute_gaussian_latitudes_and_weights.f
!
!
!     subroutine sfvpgc(nlat, nlon, isym, nt, sf, vp, idv, jdv, br, bi, cr, ci, 
!    +                   mdb, ndb, wshsgc, lshsgc, work, lwork, ierror)
!
!     given the vector spherical harmonic coefficients br, bi, cr, ci, 
!     computed by subroutine vhagc for a vector field (v, w), sfvpgc
!     computes a scalar stream function sf and scalar velocity potential
!     vp for (v, w).  (v, w) is expressed in terms of sf and vp by the
!     helmholtz relations (in mathematical spherical coordinates):
!
!          v = -1/sint*d(vp)/dlambda + d(st)/dtheta
!
!          w =  1/sint*d(st)/dlambda + d(vp)/dtheta
!
!     where sint = sin(theta).  w is the east longitudinal and v
!     is the colatitudinal component of the vector field from which
!     br, bi, cr, ci were precomputed.  required associated legendre
!     polynomials are recomputed rather than stored as they are in
!     subroutine sfvpgs. sf(i, j) and vp(i, j) are given at the i(th)
!     gaussian colatitude point theta(i) (see nlat description below)
!     and east longitude lambda(j) = (j-1)*2*pi/nlon on the sphere.
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
!            nlon = 72 for a five degree grid. nlon must be greater than
!            3.  the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!
!     isym   a parameter which determines whether the stream function and
!            velocity potential are computed on the full or half sphere
!            as follows:
!
!      = 0
!
!            the symmetries/antsymmetries described in isym=1, 2 below
!            do not exist in (v, w) about the equator.  in this case st
!            and vp are not necessarily symmetric or antisymmetric about
!            the equator.  sf and vp are computed on the entire sphere.
!            i.e., in arrays sf(i, j), vp(i, j) for i=1, ..., nlat and
!            j=1, ..., nlon.
!
!      = 1
!
!            w is antisymmetric and v is symmetric about the equator.
!            in this case sf is symmetric and vp antisymmetric about
!            the equator and are computed for the northern hemisphere
!            only.  i.e., if nlat is odd the sf(i, j), vp(i, j) are computed
!            for i=1, ..., (nlat + 1)/2 and for j=1, ..., nlon.  if nlat is
!            even then sf(i, j), vp(i, j) are computed for i=1, ..., nlat/2
!            and j=1, ..., nlon.
!
!      = 2
!
!            w is symmetric and v is antisymmetric about the equator.
!            in this case sf is antisymmetric and vp symmetric about
!            the equator and are computed for the northern hemisphere
!            only.  i.e., if nlat is odd the sf(i, j), vp(i, j) are computed
!            for i=1, ..., (nlat + 1)/2 and for j=1, ..., nlon.  if nlat is
!            even then sf(i, j), vp(i, j) are computed for i=1, ..., nlat/2
!            and j=1, ..., nlon.
!
!     nt     nt is the number of scalar and vector fields.  some
!            computational efficiency is obtained for multiple fields. arrays
!            can be three dimensional corresponding to an indexed multiple
!            vector field.  in this case multiple scalar synthesis will
!            be performed to compute sf, vp for each field.  the
!            third index is the synthesis index which assumes the values
!            k=1, ..., nt.  for a single synthesis set nt = 1.  the
!            description of the remaining parameters is simplified by
!            assuming that nt=1 or that all the arrays are two dimensional.
!
!     idv    the first dimension of the arrays sf, vp as it appears in
!            the program that calls sfvpgc. if isym = 0 then idv
!            must be at least nlat.  if isym = 1 or 2 and nlat is
!            even then idv must be at least nlat/2. if isym = 1 or 2
!            and nlat is odd then idv must be at least (nlat + 1)/2.
!
!     jdv    the second dimension of the arrays sf, vp as it appears in
!            the program that calls sfvpgc. jdv must be at least nlon.
!
!     br, bi, two or three dimensional arrays (see input parameter nt)
!     cr, ci  that contain vector spherical harmonic coefficients
!            of the vector field (v, w) as computed by subroutine vhagc.
!
!     mdb    the first dimension of the arrays br, bi, cr, ci as it
!            appears in the program that calls sfvpgc. mdb must be at
!            least min(nlat, nlon/2) if nlon is even or at least
!            min(nlat, (nlon + 1)/2) if nlon is odd.
!
!     ndb    the second dimension of the arrays br, bi, cr, ci as it
!            appears in the program that calls sfvpgc. ndb must be at
!            least nlat.
!
!     wshsgc an array which must be initialized by subroutine shsgci.
!            once initialized, wshsgc can be used repeatedly by sfvpgc
!            as long as nlon and nlat remain unchanged.  wshsgc must
!            not bel altered between calls of sfvpgc.
!
!
!     lshsgc the dimension of the array wshsgc as it appears in the
!            program that calls sfvpgc. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat + 1)/2    if nlat is odd
!
!            then lshsgc must be at least
!
!               nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls sfvpgc. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2                    if nlat is even or
!               l2 = (nlat + 1)/2                if nlat is odd
!
!            if isym is zero then lwork must be at least
!
!                nlat*((nt*nlon+max(3*l2, nlon)) + 2*l1*nt+1)
!
!            if isym is not zero then lwork must be at least
!
!                l2*(nt*nlon+max(3*nlat, nlon)) + nlat*(2*l1*nt+1)
!
!     **************************************************************
!
!     output parameters
!
!    sf, vp  two or three dimensional arrays (see input parameter nt)
!           that contains the stream function and velocity potential
!           of the vector field (v, w) whose coefficients br, bi, cr, ci
!           where precomputed by subroutine vhagc.  sf(i, j), vp(i, j)
!           are given at the i(th) gaussian colatitude point theta(i)
!           and longitude point lambda(j) = (j-1)*2*pi/nlon.  the index
!           ranges are defined above at the input parameter isym.
!
!
!    ierror = 0  no errors
!           = 1  error in the specification of nlat
!           = 2  error in the specification of nlon
!           = 3  error in the specification of isym
!           = 4  error in the specification of nt
!           = 5  error in the specification of idv
!           = 6  error in the specification of jdv
!           = 7  error in the specification of mdb
!           = 8  error in the specification of ndb
!           = 9  error in the specification of lshsgc
!           = 10 error in the specification of lwork
! **********************************************************************
!
module module_sfvpgc

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use scalar_synthesis_routines, only: &
        shsgc

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: sfvpgc

contains

    subroutine sfvpgc(nlat, nlon, isym, nt, sf, vp, idv, jdv, br, bi, cr, ci, &
        mdb, ndb, wshsgc, lshsgc, work, lwork, ierror)

        integer(ip) :: nlat, nlon, isym, nt, idv, jdv, mdb, ndb, lshsgc, lwork, ierror
        real(wp) :: sf(idv, jdv, nt), vp(idv, jdv, nt)
        real(wp) :: br(mdb, ndb, nt), bi(mdb, ndb, nt)
        real(wp) :: cr(mdb, ndb, nt), ci(mdb, ndb, nt)
        real(wp) :: wshsgc(lshsgc), work(lwork)
        integer(ip) :: imid, mmax, lzz1, labc, ls, nln, mab, mn, ia, ib, is, lwk, iwk, lwmin
        integer(ip) :: l1, l2
        !
        ! Check calling arguments
        !
        ierror = 1
        if (nlat < 3) return
        ierror = 2
        if (nlon < 4) return
        ierror = 3
        if (isym < 0 .or. isym > 2) return
        ierror = 4
        if (nt < 0) return
        ierror = 5
        imid = (nlat + 1)/2
        if ((isym == 0 .and. idv<nlat) .or. &
            (isym>0 .and. idv<imid)) return
        ierror = 6
        if (jdv < nlon) return
        ierror = 7
        mmax = min(nlat, (nlon+2)/2)
        if (mdb < min(nlat, (nlon + 1)/2)) return
        ierror = 8
        if (ndb < nlat) return
        ierror = 9
        !
        !     verify saved workspace (same as shsgc)
        !
        imid = (nlat + 1)/2
        lzz1 = 2*nlat*imid
        labc = 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
        lwmin = lzz1+labc+nlon+15
        l2 = (nlat + 1)/2
        l1 = min((nlon+2)/2, nlat)
        if (lshsgc < nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15)return
        !
        !     verify unsaved workspace (add to what shsgc requires)
        !
        ierror = 10
        ls = nlat
        if (isym> 0) ls = imid
        nln = nt*ls*nlon
        !
        !     set first dimension for a, b (as required by shsgc)
        !
        mab = min(nlat, nlon/2+1)
        mn = mab*nlat*nt
        if (lwork < nln+max(ls*nlon, 3*nlat*imid)+2*mn+nlat) return
        ierror = 0
        !
        ! Set workspace pointer indices
        !
        ia = 1
        ib = ia+mn
        is = ib+mn
        iwk = is+nlat
        lwk = lwork-2*mn-nlat
        call sfvpgc_lower_utility_routine(nlat, nlon, isym, nt, sf, vp, idv, jdv, br, bi, cr, ci, mdb, ndb, &
            work(ia), work(ib), mab, work(is), wshsgc, lshsgc, work(iwk), lwk, &
            ierror)

    contains

        subroutine sfvpgc_lower_utility_routine(nlat, nlon, isym, nt, sf, vp, idv, jdv, br, bi, cr, ci, &
            mdb, ndb, a, b, mab, fnn, wshsgc, lshsgc, wk, lwk, ierror)

            integer(ip) :: nlat, nlon, isym, nt, idv, jdv, mdb, ndb, mab, lshsgc, lwk, ierror
            real(wp) :: sf(idv, jdv, nt), vp(idv, jdv, nt)
            real(wp) :: br(mdb, ndb, nt), bi(mdb, ndb, nt), cr(mdb, ndb, nt), ci(mdb, ndb, nt)
            real(wp) :: a(mab, nlat, nt), b(mab, nlat, nt)
            real(wp) :: wshsgc(lshsgc), wk(lwk), fnn(nlat)
            integer(ip) :: n, m, mmax, k
            !
            !     set coefficient multiplyers
            !
            do n=2, nlat
                fnn(n) = 1.0/sqrt(real(n*(n-1)))
            end do
            mmax = min(nlat, (nlon + 1)/2)
            !
            !     compute sf scalar coefficients from cr, ci
            !
            do k=1, nt
                do n=1, nlat
                    do m=1, mab
                        a(m, n, k) = 0.0
                        b(m, n, k) = 0.0
                    end do
                end do
                !
                ! Compute m=0 coefficients
                !
                do n=2, nlat
                    a(1, n, k) =-fnn(n)*cr(1, n, k)
                    b(1, n, k) =-fnn(n)*ci(1, n, k)
                end do
                !
                !     compute m>0 coefficients using vector spherepack value for mmax
                !
                do m=2, mmax
                    do n=m, nlat
                        a(m, n, k) =-fnn(n)*cr(m, n, k)
                        b(m, n, k) =-fnn(n)*ci(m, n, k)
                    end do
                end do
            end do
            !
            !     synthesize a, b into st
            !
            call shsgc(nlat, nlon, isym, nt, sf, idv, jdv, a, b, mab, nlat, wshsgc, ierror)
            !
            !    set coefficients for vp from br, bi
            !
            do k=1, nt
                do n=1, nlat
                    do m=1, mab
                        a(m, n, k) = 0.0
                        b(m, n, k) = 0.0
                    end do
                end do
                 !
                 ! Compute m=0 coefficients
                 !
                do n=2, nlat
                    a(1, n, k) = fnn(n)*br(1, n, k)
                    b(1, n, k) = fnn(n)*bi(1, n, k)
                end do
                !
                !     compute m>0 coefficients using vector spherepack value for mmax
                !
                mmax = min(nlat, (nlon + 1)/2)
                do m=2, mmax
                    do n=m, nlat
                        a(m, n, k) = fnn(n)*br(m, n, k)
                        b(m, n, k) = fnn(n)*bi(m, n, k)
                    end do
                end do
            end do
            !
            !     synthesize a, b into vp
            !
            call shsgc(nlat, nlon, isym, nt, vp, idv, jdv, a, b, mab, nlat, wshsgc, ierror)

        end subroutine sfvpgc_lower_utility_routine

    end subroutine sfvpgc

end module module_sfvpgc
