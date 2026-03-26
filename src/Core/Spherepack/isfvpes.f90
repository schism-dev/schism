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
! ... file isfvpes.f
!
!     this file includes documentation and code for
!     subroutine isfvpes          i
!
! ... files which must be loaded with isfvpes.f
!
!     type_SpherepackUtility.f, type_RealPeriodicFastFourierTransform.f, vhses.f, shaes.f
!
!
!     subroutine isfvpes(nlat, nlon, isym, nt, sf, vp, idv, jdv, as, bs, av, bv, 
!    +                   mdb, ndb, wvhses, lvhses, work, lwork, ierror)
!
!     given the scalar spherical harmonic coefficients as, bs precomputed
!     by shaes for the scalar stream function sf and av, bv precomputed by
!     shaes for the scalar velocity potenital vp, subroutine isfvpes computes
!     the vector field (v, w) corresponding to sf and vp.  w is the east
!     longitudinal and v is the colatitudinal component of the vector field.
!     (v, w) is expressed in terms of sf, vp by the helmholtz relations (in
!     mathematical spherical coordinates):
!
!          v = -1/sin(theta)*d(vp)/dlambda + d(st)/dtheta
!
!          w =  1/sin(theta)*d(st)/dlambda + d(vp)/dtheta
!
!     required legendre functions are stored rather than recomputed as
!     they are in subroutine isfvpes.  v(i, j) and w(i, j) are given at
!     colatitude
!
!            theta(i) = (i-1)*pi/(nlat-1)
!
!     and east longitude
!
!            lambda(j) = (j-1)*2*pi/nlon
!
!     on the sphere.
!
!
!     input parameters
!
!     nlat   the number of colatitudes on the full sphere including the
!            poles. for example, nlat = 37 for a five degree grid.
!            nlat determines the grid increment in colatitude as
!            pi/(nlat-1).  if nlat is odd the equator is located at
!            grid point i=(nlat + 1)/2. if nlat is even the equator is
!            located half way between points i=nlat/2 and i=nlat/2+1.
!            nlat must be at least 3. note: on the half sphere, the
!            number of grid points in the colatitudinal direction is
!            nlat/2 if nlat is even or (nlat + 1)/2 if nlat is odd.
!
!     nlon   the number of distinct londitude points.  nlon determines
!            the grid increment in longitude as 2*pi/nlon. for example
!            nlon = 72 for a five degree grid. nlon must be greater than
!            3.  the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!
!     isym   a parameter which determines whether the vector field is
!            computed on the full or half sphere as follows:
!
!      = 0
!
!            the symmetries/antsymmetries described in isym=1, 2 below
!            do not exist in sf, vp about the equator.  in this case v
!            and w are not necessarily symmetric or antisymmetric about
!            equator.  v and w are computed on the entire sphere.
!            i.e., in arrays sf(i, j), vp(i, j) for i=1, ..., nlat and
!            j=1, ..., nlon.
!
!      = 1
!
!            vp is antisymmetric and sf is symmetric about the equator.
!            in this case v is symmetric and w antisymmetric about
!            the equator and are computed for the northern hemisphere
!            only.  i.e., if nlat is odd the v(i, j), w(i, j) are computed
!            for i=1, ..., (nlat + 1)/2 and for j=1, ..., nlon.  if nlat is
!            even then v(i, j), w(i, j) are computed for i=1, ..., nlat/2
!            and j=1, ..., nlon.
!
!      = 2
!
!            vp is symmetric and sf is antisymmetric about the equator.
!            in this case v is antisymmetric and w symmetric about
!            the equator and are computed for the northern hemisphere
!            only.  i.e., if nlat is odd the v(i, j), w(i, j) are computed
!            for i=1, ..., (nlat + 1)/2 and for j=1, ..., nlon.  if nlat is
!            even then v(i, j), w(i, j) are computed for i=1, ..., nlat/2
!            and j=1, ..., nlon.
!
!     nt     nt is the number of scalar and vector fields.  some
!            computational efficiency is obtained for multiple fields. arrays
!            can be three dimensional corresponding to an indexed multiple
!            vector field.  in this case multiple vector synthesis will
!            be performed to compute (v, w) for each field.  the
!            third index is the synthesis index which assumes the values
!            k=1, ..., nt.  for a single synthesis set nt = 1.  the
!            description of the remaining parameters is simplified by
!            assuming that nt=1 or that all the arrays are two dimensional.
!
!     idv    the first dimension of the arrays v, w as it appears in
!            the program that calls isfvpes. if isym = 0 then idv
!            must be at least nlat.  if isym = 1 or 2 and nlat is
!            even then idv must be at least nlat/2. if isym = 1 or 2
!            and nlat is odd then idv must be at least (nlat + 1)/2.
!
!     jdv    the second dimension of the arrays v, w as it appears in
!            the program that calls isfvpes. jdv must be at least nlon.
!
!     as, bs  two or three dimensional arrays (see input parameter nt)
!            that contain the spherical harmonic coefficients of
!            the scalar field sf as computed by subroutine shaes.
!
!     av, bv  two or three dimensional arrays (see input parameter nt)
!            that contain the spherical harmonic coefficients of
!            the scalar field vp as computed by subroutine shaes.
!
!     mdb    the first dimension of the arrays as, bs, av, bv as it
!            appears in the program that calls isfvpes. mdb must be at
!            least min(nlat, (nlon+2)/2) if nlon is even or at least
!            min(nlat, (nlon + 1)/2) if nlon is odd.
!
!     ndb    the second dimension of the arrays as, bs, av, bv as it
!            appears in the program that calls isfvpes. ndb must be at
!            least nlat.
!
!     wvhses an array which must be initialized by subroutine vhsesi.
!            once initialized, wvhses can be used repeatedly by isfvpes
!            as long as nlon and nlat remain unchanged.  wvhses must
!            not bel altered between calls of isfvpes.
!
!
!     lvhses the dimension of the array wvhses as it appears in the
!            program that calls isfvpes. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat + 1)/2    if nlat is odd
!
!            then lvhses must be at least
!
!               l1*l2*(2*nlat-l1+1)+nlon+15
!
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls isfvpes. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2                    if nlat is even or
!               l2 = (nlat + 1)/2                if nlat is odd
!
!            if isym=0 then lwork must be at least
!
!                       nlat*((2*nt+1)*nlon + 4*l1*nt + 1)
!
!            if isym=1 or 2 then lwork must be at least
!
!                       (2*nt+1)*l2*nlon + nlat*(4*l1*nt + 1)
!
!
!     **************************************************************
!
!     output parameters
!
!    v, w    two or three dimensional arrays (see input parameter nt)
!           that contains the vector field corresponding to the stream
!           function sf and velocity potential vp whose coefficients, 
!           as, bs (for sf) and av, bv (for vp), were precomputed by
!           subroutine shaes.  v(i, j) and w(i, j) are given at the
!           colatitude point
!
!                theta(i) = (i-1)*pi/(nlat-1)
!
!           and longitude point
!
!                lambda(j) = (j-1)*2*pi/nlon
!
!           the index ranges are defined above at the input parameter isym.
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
!           = 9  error in the specification of lvhses
!           = 10 error in the specification of lwork
! **********************************************************************
!
module module_isfvpes

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use vector_synthesis_routines, only: &
        vhses

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: isfvpes

contains

    subroutine isfvpes(nlat, nlon, isym, nt, v, w, idv, jdv, as, bs, av, bv, &
        mdb, ndb, wvhses, lvhses, work, lwork, ierror)

        integer(ip) :: nlat, nlon, isym, nt, idv, jdv, mdb, ndb, lvhses, lwork, ierror
        real(wp) :: v(idv, jdv, nt), w(idv, jdv, nt)
        real(wp) :: as(mdb, ndb, nt), bs(mdb, ndb, nt)
        real(wp) :: av(mdb, ndb, nt), bv(mdb, ndb, nt)
        real(wp) :: wvhses(lvhses), work(lwork)
        integer(ip) :: l1, l2, mn, is, lwk, iwk
        integer(ip) :: ibr, ibi, icr, ici
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
        l2 = (nlat + 1)/2
        if ((isym == 0 .and. idv<nlat) .or. &
            (isym>0 .and. idv<l2)) return
        ierror = 6
        if (jdv < nlon) return
        ierror = 7
        l1 = min(nlat, (nlon + 1)/2)
        if (mdb < min(nlat, (nlon+2)/2)) return
        ierror = 8
        if (ndb < nlat) return
        ierror = 9
        l1 = min(nlat, (nlon+2)/2)
        if (lvhses < l1*l2*(2*nlat-l1+1)+nlon+15) return
        ierror = 10
        if (isym == 0) then
            if (lwork < nlat*((2*nt+1)*nlon+4*l1*nt+1)) return
        else
            if (lwork < (2*nt+1)*nlon+nlat*(4*l1*nt+1)) return
        end if
        !
        !     set first dimension for br, bi, cr, ci (as required by vhses)
        !
        mn = l1*nlat*nt
        ierror = 0
        !
        ! Set workspace pointer indices
        !
        ibr = 1
        ibi = ibr+mn
        icr = ibi+mn
        ici = icr+mn
        is = ici+mn
        iwk = is+nlat
        lwk = lwork-4*mn-nlat
        call isfvpes1(nlat, nlon, isym, nt, v, w, idv, jdv, as, bs, av, bv, mdb, &
            ndb, work(ibr), work(ibi), work(icr), work(ici), l1, work(is), &
            wvhses, lvhses, work(iwk), lwk, ierror)

    contains

        subroutine isfvpes1(nlat, nlon, isym, nt, v, w, idv, jdv, as, bs, av, bv, &
            mdb, ndb, br, bi, cr, ci, mab, fnn, wvhses, lvhses, wk, lwk, ierror)

            integer(ip) :: nlat, nlon, isym, nt, idv, jdv, mdb, ndb, mab, lvhses, lwk, ierror
            real(wp) :: v(idv, jdv, nt), w(idv, jdv, nt)
            real(wp) :: as(mdb, ndb, nt), bs(mdb, ndb, nt)
            real(wp) :: av(mdb, ndb, nt), bv(mdb, ndb, nt)
            real(wp) :: br(mab, nlat, nt), bi(mab, nlat, nt)
            real(wp) :: cr(mab, nlat, nt), ci(mab, nlat, nt)
            real(wp) :: wvhses(lvhses), wk(lwk), fnn(nlat)
            integer(ip) :: n, m, mmax, k, ityp
            !
            !     set coefficient multiplyers
            !
            do n=2, nlat
                fnn(n) = -sqrt(real(n*(n-1)))
            end do
            mmax = min(nlat, (nlon + 1)/2)
            !
            !     compute (v, w) coefficients from as, bs, av, bv
            !
            do k=1, nt
                do n=1, nlat
                    do m=1, mab
                        br(m, n, k) = 0.0
                        bi(m, n, k) = 0.0
                        cr(m, n, k) = 0.0
                        ci(m, n, k) = 0.0
                    end do
                end do
                !
                ! Compute m=0 coefficients
                !
                do n=2, nlat
                    br(1, n, k) = -fnn(n)*av(1, n, k)
                    bi(1, n, k) = -fnn(n)*bv(1, n, k)
                    cr(1, n, k) =  fnn(n)*as(1, n, k)
                    ci(1, n, k) =  fnn(n)*bs(1, n, k)
                end do
                !
                !     compute m>0 coefficients using vector spherepack value for mmax
                !
                do m=2, mmax
                    do n=m, nlat
                        br(m, n, k) = -fnn(n)*av(m, n, k)
                        bi(m, n, k) = -fnn(n)*bv(m, n, k)
                        cr(m, n, k) =  fnn(n)*as(m, n, k)
                        ci(m, n, k) =  fnn(n)*bs(m, n, k)
                    end do
                end do
            end do
            !
            !     synthesize br, bi, cr, ci into (v, w)
            !
            select case (isym)
                case (0)
                    ityp = 0
                case (1)
                    ityp = 3
                case (2)
                    ityp = 6
            end select

            call vhses(nlat, nlon, ityp, nt, v, w, idv, jdv, br, bi, cr, ci, &
                mab, nlat, wvhses, ierror)

        end subroutine isfvpes1

    end subroutine isfvpes

end module module_isfvpes
