!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                         Spherepack                            *
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
submodule(colatitudinal_derivative_routines) colatitudinal_derivative_regular_grid_saved

contains

    !     subroutine vtses(nlat, nlon, ityp, nt, vt, wt, idvw, jdvw, br, bi, cr, ci, &
    !                      mdab, ndab, wvts, lwvts, work, lwork, ierror)
    !
    !     given the vector harmonic analysis br, bi, cr, and ci (computed
    !     by subroutine vhaes) of some vector function (v, w), this
    !     subroutine computes the vector function (vt, wt) which is
    !     the derivative of (v, w) with respect to colatitude theta. vtses
    !     is similar to vhses except the vector harmonics are replaced by
    !     their derivative with respect to colatitude with the result that
    !     (vt, wt) is computed instead of (v, w). vt(i, j) is the derivative
    !     of the colatitudinal component v(i, j) at the point theta(i) =
    !     (i-1)*pi/(nlat-1) and longitude phi(j) = (j-1)*2*pi/nlon. the
    !     spectral representation of (vt, wt) is given below at output
    !     parameters vt, wt.
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
    !            nlat/2 if nlat is even or(nlat + 1)/2 if nlat is odd.
    !
    !     nlon   the number of distinct londitude points.  nlon determines
    !            the grid increment in longitude as 2*pi/nlon. for example
    !            nlon = 72 for a five degree grid. nlon must be greater
    !            than zero. the axisymmetric case corresponds to nlon=1.
    !            the efficiency of the computation is improved when nlon
    !            is a product of small prime numbers.
    !
    !     ityp   = 0  no symmetries exist about the equator. the synthesis
    !                 is performed on the entire sphere. i.e. the arrays
    !                 vt(i, j), wt(i, j) are computed for i=1, ..., nlat and
    !                 j=1, ..., nlon.
    !
    !            = 1  no symmetries exist about the equator however the
    !                 the coefficients cr and ci are zero. the synthesis
    !                 is performed on the entire sphere.  i.e. the arrays
    !                 vt(i, j), wt(i, j) are computed for i=1, ..., nlat and
    !                 j=1, ..., nlon.
    !
    !            = 2  no symmetries exist about the equator however the
    !                 the coefficients br and bi are zero. the synthesis
    !                 is performed on the entire sphere.  i.e. the arrays
    !                 vt(i, j), wt(i, j) are computed for i=1, ..., nlat and
    !                 j=1, ..., nlon.
    !
    !            = 3  vt is odd and wt is even about the equator. the
    !                 synthesis is performed on the northern hemisphere
    !                 only.  i.e., if nlat is odd the arrays vt(i, j), wt(i, j)
    !                 are computed for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon.
    !                 if nlat is even the arrays vt(i, j), wt(i, j) are computed
    !                 for i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !            = 4  vt is odd and wt is even about the equator and the
    !                 coefficients cr and ci are zero. the synthesis is
    !                 performed on the northern hemisphere only. i.e. if
    !                 nlat is odd the arrays vt(i, j), wt(i, j) are computed
    !                 for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon. if nlat is
    !                 even the arrays vt(i, j), wt(i, j) are computed for
    !                 i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !            = 5  vt is odd and wt is even about the equator and the
    !                 coefficients br and bi are zero. the synthesis is
    !                 performed on the northern hemisphere only. i.e. if
    !                 nlat is odd the arrays vt(i, j), wt(i, j) are computed
    !                 for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon. if nlat is
    !                 even the arrays vt(i, j), wt(i, j) are computed for
    !                 i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !            = 6  vt is even and wt is odd about the equator. the
    !                 synthesis is performed on the northern hemisphere
    !                 only.  i.e., if nlat is odd the arrays vt(i, j), wt(i, j)
    !                 are computed for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon.
    !                 if nlat is even the arrays vt(i, j), wt(i, j) are computed
    !                 for i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !            = 7  vt is even and wt is odd about the equator and the
    !                 coefficients cr and ci are zero. the synthesis is
    !                 performed on the northern hemisphere only. i.e. if
    !                 nlat is odd the arrays vt(i, j), wt(i, j) are computed
    !                 for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon. if nlat is
    !                 even the arrays vt(i, j), wt(i, j) are computed for
    !                 i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !            = 8  vt is even and wt is odd about the equator and the
    !                 coefficients br and bi are zero. the synthesis is
    !                 performed on the northern hemisphere only. i.e. if
    !                 nlat is odd the arrays vt(i, j), wt(i, j) are computed
    !                 for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon. if nlat is
    !                 even the arrays vt(i, j), wt(i, j) are computed for
    !                 i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !     nt     the number of syntheses.  in the program that calls vtses,
    !            the arrays vt, wt, br, bi, cr, and ci can be three dimensional
    !            in which case multiple syntheses will be performed.
    !            the third index is the synthesis index which assumes the
    !            values k=1, ..., nt.  for a single synthesis set nt=1. the
    !            discription of the remaining parameters is simplified
    !            by assuming that nt=1 or that all the arrays are two
    !            dimensional.
    !
    !     idvw   the first dimension of the arrays vt, wt as it appears in
    !            the program that calls vtses. if ityp <= 2 then idvw
    !            must be at least nlat.  if ityp > 2 and nlat is
    !            even then idvw must be at least nlat/2. if ityp > 2
    !            and nlat is odd then idvw must be at least(nlat + 1)/2.
    !
    !     jdvw   the second dimension of the arrays vt, wt as it appears in
    !            the program that calls vtses. jdvw must be at least nlon.
    !
    !     br, bi  two or three dimensional arrays (see input parameter nt)
    !     cr, ci  that contain the vector spherical harmonic coefficients
    !            of (v, w) as computed by subroutine vhaes.
    !
    !     mdab   the first dimension of the arrays br, bi, cr, and ci as it
    !            appears in the program that calls vtses. mdab must be at
    !            least min(nlat, nlon/2) if nlon is even or at least
    !            min(nlat, (nlon + 1)/2) if nlon is odd.
    !
    !     ndab   the second dimension of the arrays br, bi, cr, and ci as it
    !            appears in the program that calls vtses. ndab must be at
    !            least nlat.
    !
    !     wvts   an array which must be initialized by subroutine vtsesi.
    !            once initialized, wvts can be used repeatedly by vtses
    !            as long as nlon and nlat remain unchanged.  wvts must
    !            not be altered between calls of vtses.
    !
    !     lwvts  the dimension of the array wvts as it appears in the
    !            program that calls vtses. define
    !
    !               l1 = min(nlat, nlon/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 =(nlat + 1)/2    if nlat is odd
    !
    !            then lwvts must be at least
    !
    !                 l1*l2*(2*nlat-l1+1)+nlon+15
    !
    !     output parameters
    !
    !     vt, wt  two or three dimensional arrays (see input parameter nt)
    !            in which the derivative of (v, w) with respect to
    !            colatitude theta is stored. vt(i, j), wt(i, j) contain the
    !            derivatives at colatitude theta(i) = (i-1)*pi/(nlat-1)
    !            and longitude phi(j) = (j-1)*2*pi/nlon. the index ranges
    !            are defined above at the input parameter ityp. vt and wt
    !            are computed from the formulas for v and w given in
    !            subroutine vhses but with vbar and wbar replaced with
    !            their derivatives with respect to colatitude. these
    !            derivatives are denoted by vtbar and wtbar.
    !
    !
    !   *************************************************************
    !
    !   in terms of real variables this expansion takes the form
    !
    !             for i=1, ..., nlat and  j=1, ..., nlon
    !
    !     vt(i, j) = the sum from n=1 to n=nlat-1 of
    !
    !               .5*br(1, n+1)*vtbar(0, n, theta(i))
    !
    !     plus the sum from m=1 to m=mmax-1 of the sum from n=m to
    !     n=nlat-1 of the real part of
    !
    !       (br(m+1, n+1)*vtbar(m, n, theta(i))
    !                   -ci(m+1, n+1)*wtbar(m, n, theta(i)))*cos(m*phi(j))
    !      -(bi(m+1, n+1)*vtbar(m, n, theta(i))
    !                   +cr(m+1, n+1)*wtbar(m, n, theta(i)))*sin(m*phi(j))
    !
    !    and for i=1, ..., nlat and  j=1, ..., nlon
    !
    !     wt(i, j) = the sum from n=1 to n=nlat-1 of
    !
    !              -.5*cr(1, n+1)*vtbar(0, n, theta(i))
    !
    !     plus the sum from m=1 to m=mmax-1 of the sum from n=m to
    !     n=nlat-1 of the real part of
    !
    !      -(cr(m+1, n+1)*vtbar(m, n, theta(i))
    !                   +bi(m+1, n+1)*wtbar(m, n, theta(i)))*cos(m*phi(j))
    !      +(ci(m+1, n+1)*vtbar(m, n, theta(i))
    !                   -br(m+1, n+1)*wtbar(m, n, theta(i)))*sin(m*phi(j))
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
    !            = 9  error in the specification of lwvts
    !
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


        ! Local variables
        integer(ip) :: idv
        integer(ip) :: idz
        integer(ip) :: imid
        integer(ip) :: ist
        integer(ip) :: iw1
        integer(ip) :: iw2
        integer(ip) :: iw3
        integer(ip) :: iw4
        integer(ip) :: jw1
        integer(ip) :: jw2
        integer(ip) :: lnl
        integer(ip) :: lzimn
        integer(ip) :: mmax, lwork, lwvts

        lwvts = size(wvts)

        ! Check calling arguments
        ierror = 1
        if(nlat < 3) return
        ierror = 2
        if (nlon < 1) return
        ierror = 3
        if (ityp<0 .or. ityp>8) return
        ierror = 4
        if (nt < 0) return
        ierror = 5
        imid =(nlat + 1)/2
        if ((ityp<=2 .and. idvw<nlat) .or. &
            (ityp>2 .and. idvw<imid)) return
        ierror = 6
        if (jdvw < nlon) return
        ierror = 7
        mmax = min(nlat, (nlon + 1)/2)
        if (mdab < mmax) return
        ierror = 8
        if (ndab < nlat) return
        ierror = 9
        idz = (mmax*(2*nlat-mmax+1))/2
        lzimn = idz*imid
        if (lwvts < 2*lzimn+nlon+15) return
        ierror = 10
        idv = nlat
        if (ityp > 2) idv = imid
        lnl = nt*idv*nlon
        lwork = 2*lnl+idv*nlon
        ierror = 0
        ist = 0
        if (ityp <= 2) ist = imid

        block
            real(wp) :: work(lwork)
            iw1 = ist+1
            iw2 = lnl+1
            iw3 = iw2+ist
            iw4 = iw2+lnl
            jw1 = lzimn+1
            jw2 = jw1+lzimn

            call vtses_lower_utility_routine(nlat, nlon, ityp, nt, imid, idvw, jdvw, vt, wt, mdab, ndab, &
                br, bi, cr, ci, idv, work, work(iw1:), work(iw2:), work(iw3:), &
                work(iw4:), idz, wvts, wvts(jw1:), wvts(jw2:))
        end block

    end subroutine vtses

    !     subroutine vtsesi(nlat, nlon, wvts, ierror)
    !
    !     subroutine vtsesi initializes the array wvts which can then be
    !     used repeatedly by subroutine vtses until nlat or nlon is changed.
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
    !            nlat/2 if nlat is even or(nlat + 1)/2 if nlat is odd.
    !
    !     nlon   the number of distinct londitude points.  nlon determines
    !            the grid increment in longitude as 2*pi/nlon. for example
    !            nlon = 72 for a five degree grid. nlon must be greater
    !            than zero. the axisymmetric case corresponds to nlon=1.
    !            the efficiency of the computation is improved when nlon
    !            is a product of small prime numbers.
    !
    !     lwvts  the dimension of the array wvts as it appears in the
    !            program that calls vtses. define
    !
    !               l1 = min(nlat, nlon/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 =(nlat + 1)/2    if nlat is odd
    !
    !            then lwvts must be at least
    !
    !                  l1*l2*(2*nlat-l1+1)+nlon+15
    !
    !     output parameters
    !
    !     wvts   an array which is initialized for use by subroutine vtses.
    !            once initialized, wvts can be used repeatedly by vtses
    !            as long as nlat or nlon remain unchanged.  wvts must not
    !            be altered between calls of vtses.
    !
    !
    !     ierror = 0  no errors
    !            = 1  error in the specification of nlat
    !            = 2  error in the specification of nlon
    !            = 3  error in the specification of lwvts
    !
    module subroutine vtsesi(nlat, nlon, wvts, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(out) :: wvts(:)
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip)             :: imid, labc, lzimn
        integer(ip)             :: mmax, lwork, ldwork
        type(SpherepackUtility) :: util

        ! Check calling arguments
        call check_init_calling_arguments(&
            nlat, nlon, wvts, ierror, get_lwvts_saved(nlat,nlon))

        ! Check error flag
        if (ierror /= 0) return

        mmax = min(nlat, (nlon + 1)/2)
        imid =(nlat + 1)/2
        lzimn = (imid*mmax*(2*nlat-mmax+1))/2
        labc = 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2

        ! Set required workspace size
        lwork = 5*nlat*imid+labc
        ldwork = 2*(nlat + 1)

        block
            real(wp)    :: work(lwork), dwork(ldwork)
            integer(ip) :: iw1, iw2, iw3, idz

            ! Set worspace index pointers
            iw1 = 3*nlat*imid+1
            iw2 = lzimn+1
            iw3 = 2*lzimn+1
            idz = (mmax*(2*nlat-mmax+1))/2

            call vtsesi_lower_utility_routine(nlat, nlon, imid, wvts, wvts(iw2:), &
                idz, work, work(iw1:), dwork)

            call util%hfft%initialize(nlon, wvts(iw3:))
        end block

    end subroutine vtsesi

    subroutine vtses_lower_utility_routine(nlat, nlon, ityp, nt, imid, idvw, jdvw, vt, wt, mdab, &
        ndab, br, bi, cr, ci, idv, vte, vto, wte, wto, work, idz, vb, wb, wrfft)

        real(wp) :: bi
        real(wp) :: br
        real(wp) :: ci
        real(wp) :: cr
        integer(ip) :: idv
        integer(ip) :: idvw
        integer(ip) :: idz
        integer(ip) :: imid
        integer(ip) :: imm1
        integer(ip) :: ityp
        integer(ip) :: jdvw
        integer(ip) :: m
        integer(ip) :: mb
        integer(ip) :: mdab
        integer(ip) :: mlat
        integer(ip) :: mlon
        integer(ip) :: mmax
        integer(ip) :: mn
        integer(ip) :: mp1
        integer(ip) :: mp2
        integer(ip) :: ndab
        integer(ip) :: ndo1
        integer(ip) :: ndo2
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: nlp1
        integer(ip) :: np1
        integer(ip) :: nt
        real(wp) :: vb
        real(wp) :: vt
        real(wp) :: vte
        real(wp) :: vto
        real(wp) :: wb
        real(wp) :: work
        real(wp) :: wrfft
        real(wp) :: wt
        real(wp) :: wte
        real(wp) :: wto
        dimension vt(idvw, jdvw, nt), wt(idvw, jdvw, nt), br(mdab, ndab, nt), &
            bi(mdab, ndab, nt), cr(mdab, ndab, nt), ci(mdab, ndab, nt), &
            vte(idv, nlon, nt), vto(idv, nlon, nt), wte(idv, nlon, nt), &
            wto(idv, nlon, nt), work(*), wrfft(*), &
            vb(imid, *), wb(imid, *)

        integer(ip) :: i, j, k ! Counters

        type(SpherepackUtility) :: util

        nlp1 = nlat+1
        mlat = mod(nlat, 2)
        mlon = mod(nlon, 2)
        mmax = min(nlat, (nlon + 1)/2)

        select case(mlat)
            case(0)
                imm1 = imid
                ndo1 = nlat
                ndo2 = nlat-1
            case default
                imm1 = imid-1
                ndo1 = nlat-1
                ndo2 = nlat
        end select

        ! Preset even fields to zero
        vte = ZERO
        wte = ZERO

        vector_symmetry_cases: select case (ityp)
            case (0)
                        ! case ityp=0   no symmetries
                        !
                        ! case m = 0
                do k=1, nt
                    do np1=2, ndo2, 2
                        do i=1, imm1
                            vto(i, 1, k)=vto(i, 1, k)+br(1, np1, k)*vb(i, np1)
                            wto(i, 1, k)=wto(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                do k=1, nt
                    do np1=3, ndo1, 2
                        do i=1, imid
                            vte(i, 1, k)=vte(i, 1, k)+br(1, np1, k)*vb(i, np1)
                            wte(i, 1, k)=wte(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) exit vector_symmetry_cases
                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1
                    if (mp1 <= ndo1) then
                        do k=1, nt
                            do np1=mp1, ndo1, 2
                                mn = mb+np1
                                do i=1, imm1
                                    vte(i, 2*mp1-2, k) = vte(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                    vto(i, 2*mp1-2, k) = vto(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                    vte(i, 2*mp1-1, k) = vte(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                    vto(i, 2*mp1-1, k) = vto(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                    wte(i, 2*mp1-2, k) = wte(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                    wto(i, 2*mp1-2, k) = wto(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                    wte(i, 2*mp1-1, k) = wte(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                    wto(i, 2*mp1-1, k) = wto(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                end do
                                if (mlat /= 0) then
                                    vte(imid, 2*mp1-2, k) = vte(imid, 2*mp1-2, k) &
                                        +br(mp1, np1, k)*vb(imid, mn)
                                    vte(imid, 2*mp1-1, k) = vte(imid, 2*mp1-1, k) &
                                        +bi(mp1, np1, k)*vb(imid, mn)
                                    wte(imid, 2*mp1-2, k) = wte(imid, 2*mp1-2, k) &
                                        -cr(mp1, np1, k)*vb(imid, mn)
                                    wte(imid, 2*mp1-1, k) = wte(imid, 2*mp1-1, k) &
                                        -ci(mp1, np1, k)*vb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                    if (mp2 <= ndo2) then
                        do k=1, nt
                            do np1=mp2, ndo2, 2
                                mn = mb+np1
                                do i=1, imm1
                                    vto(i, 2*mp1-2, k) = vto(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                    vte(i, 2*mp1-2, k) = vte(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                    vto(i, 2*mp1-1, k) = vto(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                    vte(i, 2*mp1-1, k) = vte(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                    wto(i, 2*mp1-2, k) = wto(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                    wte(i, 2*mp1-2, k) = wte(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                    wto(i, 2*mp1-1, k) = wto(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                    wte(i, 2*mp1-1, k) = wte(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                end do
                                if (mlat /= 0) then
                                    vte(imid, 2*mp1-2, k) = vte(imid, 2*mp1-2, k) &
                                        -ci(mp1, np1, k)*wb(imid, mn)
                                    vte(imid, 2*mp1-1, k) = vte(imid, 2*mp1-1, k) &
                                        +cr(mp1, np1, k)*wb(imid, mn)
                                    wte(imid, 2*mp1-2, k) = wte(imid, 2*mp1-2, k) &
                                        -bi(mp1, np1, k)*wb(imid, mn)
                                    wte(imid, 2*mp1-1, k) = wte(imid, 2*mp1-1, k) &
                                        +br(mp1, np1, k)*wb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                end do
            case(1)
                        ! case ityp=1   no symmetries,  cr and ci equal zero
                        !
                        ! case m = 0
                do k=1, nt
                    do np1=2, ndo2, 2
                        do i=1, imm1
                            vto(i, 1, k)=vto(i, 1, k)+br(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                do k=1, nt
                    do np1=3, ndo1, 2
                        do i=1, imid
                            vte(i, 1, k)=vte(i, 1, k)+br(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) exit vector_symmetry_cases
                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1
                    if (mp1 <= ndo1) then
                        do k=1, nt
                            do np1=mp1, ndo1, 2
                                mn = mb+np1
                                do i=1, imm1
                                    vte(i, 2*mp1-2, k) = vte(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                    vte(i, 2*mp1-1, k) = vte(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                    wto(i, 2*mp1-2, k) = wto(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                    wto(i, 2*mp1-1, k) = wto(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                end do
                                if (mlat /= 0) then
                                    vte(imid, 2*mp1-2, k) = vte(imid, 2*mp1-2, k) &
                                        +br(mp1, np1, k)*vb(imid, mn)
                                    vte(imid, 2*mp1-1, k) = vte(imid, 2*mp1-1, k) &
                                        +bi(mp1, np1, k)*vb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                    if (mp2 <= ndo2) then
                        do k=1, nt
                            do np1=mp2, ndo2, 2
                                mn = mb+np1
                                do i=1, imm1
                                    vto(i, 2*mp1-2, k) = vto(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                    vto(i, 2*mp1-1, k) = vto(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                    wte(i, 2*mp1-2, k) = wte(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                    wte(i, 2*mp1-1, k) = wte(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                end do
                                if (mlat /= 0) then
                                    wte(imid, 2*mp1-2, k) = wte(imid, 2*mp1-2, k) &
                                        -bi(mp1, np1, k)*wb(imid, mn)
                                    wte(imid, 2*mp1-1, k) = wte(imid, 2*mp1-1, k) &
                                        +br(mp1, np1, k)*wb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                end do
            case(2)
                        ! case ityp=2   no symmetries,  br and bi are equal to zero
                        !
                        ! case m = 0
                do k=1, nt
                    do np1=2, ndo2, 2
                        do i=1, imm1
                            wto(i, 1, k)=wto(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                do k=1, nt
                    do np1=3, ndo1, 2
                        do i=1, imid
                            wte(i, 1, k)=wte(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) exit vector_symmetry_cases
                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1
                    if (mp1 <= ndo1) then
                        do k=1, nt
                            do np1=mp1, ndo1, 2
                                mn = mb+np1
                                do i=1, imm1
                                    vto(i, 2*mp1-2, k) = vto(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                    vto(i, 2*mp1-1, k) = vto(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                    wte(i, 2*mp1-2, k) = wte(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                    wte(i, 2*mp1-1, k) = wte(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                end do
                                if (mlat /= 0) then
                                    wte(imid, 2*mp1-2, k) = wte(imid, 2*mp1-2, k) &
                                        -cr(mp1, np1, k)*vb(imid, mn)
                                    wte(imid, 2*mp1-1, k) = wte(imid, 2*mp1-1, k) &
                                        -ci(mp1, np1, k)*vb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                    if (mp2 <= ndo2) then
                        do k=1, nt
                            do np1=mp2, ndo2, 2
                                mn = mb+np1
                                do i=1, imm1
                                    vte(i, 2*mp1-2, k) = vte(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                    vte(i, 2*mp1-1, k) = vte(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                    wto(i, 2*mp1-2, k) = wto(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                    wto(i, 2*mp1-1, k) = wto(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                end do
                                if (mlat /= 0) then
                                    vte(imid, 2*mp1-2, k) = vte(imid, 2*mp1-2, k) &
                                        -ci(mp1, np1, k)*wb(imid, mn)
                                    vte(imid, 2*mp1-1, k) = vte(imid, 2*mp1-1, k) &
                                        +cr(mp1, np1, k)*wb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                end do
            case(3)
                        ! case ityp=3   v odd,  w even
                        !
                do k=1, nt
                    do np1=2, ndo2, 2
                        do i=1, imm1
                            vto(i, 1, k)=vto(i, 1, k)+br(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                do k=1, nt
                    do np1=3, ndo1, 2
                        do i=1, imid
                            wte(i, 1, k)=wte(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) exit vector_symmetry_cases
                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1
                    if (mp1 <= ndo1) then
                        do k=1, nt
                            do np1=mp1, ndo1, 2
                                mn = mb+np1
                                do i=1, imm1
                                    vto(i, 2*mp1-2, k) = vto(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                    vto(i, 2*mp1-1, k) = vto(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                    wte(i, 2*mp1-2, k) = wte(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                    wte(i, 2*mp1-1, k) = wte(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                end do
                                if (mlat /= 0) then
                                    wte(imid, 2*mp1-2, k) = wte(imid, 2*mp1-2, k) &
                                        -cr(mp1, np1, k)*vb(imid, mn)
                                    wte(imid, 2*mp1-1, k) = wte(imid, 2*mp1-1, k) &
                                        -ci(mp1, np1, k)*vb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                    if (mp2 <= ndo2) then
                        do k=1, nt
                            do np1=mp2, ndo2, 2
                                mn = mb+np1
                                do i=1, imm1
                                    vto(i, 2*mp1-2, k) = vto(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                    vto(i, 2*mp1-1, k) = vto(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                    wte(i, 2*mp1-2, k) = wte(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                    wte(i, 2*mp1-1, k) = wte(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                end do
                                if (mlat /= 0) then
                                    wte(imid, 2*mp1-2, k) = wte(imid, 2*mp1-2, k) &
                                        -bi(mp1, np1, k)*wb(imid, mn)
                                    wte(imid, 2*mp1-1, k) = wte(imid, 2*mp1-1, k) &
                                        +br(mp1, np1, k)*wb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                end do
            case(4)
                        ! case ityp=4   v odd,  w even, and both cr and ci equal zero
                        !
                        ! case m = 0
                do k=1, nt
                    do np1=2, ndo2, 2
                        do i=1, imm1
                            vto(i, 1, k)=vto(i, 1, k)+br(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) exit vector_symmetry_cases
                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1
                    if (mp2 <= ndo2) then
                        do k=1, nt
                            do np1=mp2, ndo2, 2
                                mn = mb+np1
                                do i=1, imm1
                                    vto(i, 2*mp1-2, k) = vto(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                    vto(i, 2*mp1-1, k) = vto(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                    wte(i, 2*mp1-2, k) = wte(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                    wte(i, 2*mp1-1, k) = wte(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                end do
                                if (mlat /= 0) then
                                    wte(imid, 2*mp1-2, k) = wte(imid, 2*mp1-2, k) &
                                        -bi(mp1, np1, k)*wb(imid, mn)
                                    wte(imid, 2*mp1-1, k) = wte(imid, 2*mp1-1, k) &
                                        +br(mp1, np1, k)*wb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                end do
            case(5)
                        ! case ityp=5   v odd,  w even,     br and bi equal zero
                        !
                        ! case m = 0
                do k=1, nt
                    do np1=3, ndo1, 2
                        do i=1, imid
                            wte(i, 1, k)=wte(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) exit vector_symmetry_cases
                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1
                    if (mp1 <= ndo1) then
                        do k=1, nt
                            do np1=mp1, ndo1, 2
                                mn = mb+np1
                                do i=1, imm1
                                    vto(i, 2*mp1-2, k) = vto(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                    vto(i, 2*mp1-1, k) = vto(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                    wte(i, 2*mp1-2, k) = wte(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                    wte(i, 2*mp1-1, k) = wte(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                end do
                                if (mlat /= 0) then
                                    wte(imid, 2*mp1-2, k) = wte(imid, 2*mp1-2, k) &
                                        -cr(mp1, np1, k)*vb(imid, mn)
                                    wte(imid, 2*mp1-1, k) = wte(imid, 2*mp1-1, k) &
                                        -ci(mp1, np1, k)*vb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                end do
            case(6)
                        ! case ityp=6   v even  ,  w odd
                        !
                        ! case m = 0
                do k=1, nt
                    do np1=2, ndo2, 2
                        do i=1, imm1
                            wto(i, 1, k)=wto(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                do k=1, nt
                    do np1=3, ndo1, 2
                        do i=1, imid
                            vte(i, 1, k)=vte(i, 1, k)+br(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) exit vector_symmetry_cases
                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1
                    if (mp1 <= ndo1) then
                        do k=1, nt
                            do np1=mp1, ndo1, 2
                                mn = mb+np1
                                do i=1, imm1
                                    vte(i, 2*mp1-2, k) = vte(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                    vte(i, 2*mp1-1, k) = vte(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                    wto(i, 2*mp1-2, k) = wto(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                    wto(i, 2*mp1-1, k) = wto(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                end do
                                if (mlat /= 0) then
                                    vte(imid, 2*mp1-2, k) = vte(imid, 2*mp1-2, k) &
                                        +br(mp1, np1, k)*vb(imid, mn)
                                    vte(imid, 2*mp1-1, k) = vte(imid, 2*mp1-1, k) &
                                        +bi(mp1, np1, k)*vb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                    if (mp2 <= ndo2) then
                        do k=1, nt
                            do np1=mp2, ndo2, 2
                                mn = mb+np1
                                do i=1, imm1
                                    vte(i, 2*mp1-2, k) = vte(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                    vte(i, 2*mp1-1, k) = vte(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                    wto(i, 2*mp1-2, k) = wto(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                    wto(i, 2*mp1-1, k) = wto(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                end do
                                if (mlat /= 0) then
                                    vte(imid, 2*mp1-2, k) = vte(imid, 2*mp1-2, k) &
                                        -ci(mp1, np1, k)*wb(imid, mn)
                                    vte(imid, 2*mp1-1, k) = vte(imid, 2*mp1-1, k) &
                                        +cr(mp1, np1, k)*wb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                end do
            case(7)
                        ! case ityp=7   v even, w odd   cr and ci equal zero
                        !
                        ! case m = 0
                do k=1, nt
                    do np1=3, ndo1, 2
                        do i=1, imid
                            vte(i, 1, k)=vte(i, 1, k)+br(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) exit vector_symmetry_cases
                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1
                    if (mp1 <= ndo1) then
                        do k=1, nt
                            do np1=mp1, ndo1, 2
                                mn = mb+np1
                                do i=1, imm1
                                    vte(i, 2*mp1-2, k) = vte(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                    vte(i, 2*mp1-1, k) = vte(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                    wto(i, 2*mp1-2, k) = wto(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                    wto(i, 2*mp1-1, k) = wto(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                end do
                                if (mlat /= 0) then
                                    vte(imid, 2*mp1-2, k) = vte(imid, 2*mp1-2, k) &
                                        +br(mp1, np1, k)*vb(imid, mn)
                                    vte(imid, 2*mp1-1, k) = vte(imid, 2*mp1-1, k) &
                                        +bi(mp1, np1, k)*vb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                end do
            case(8)
                        ! case ityp=8   v even,  w odd,   br and bi equal zero
                        !
                        ! case m = 0
                do k=1, nt
                    do np1=2, ndo2, 2
                        do i=1, imm1
                            wto(i, 1, k)=wto(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) exit vector_symmetry_cases
                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1
                    if (mp2 <= ndo2) then
                        do k=1, nt
                            do np1=mp2, ndo2, 2
                                mn = mb+np1
                                do i=1, imm1
                                    vte(i, 2*mp1-2, k) = vte(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                    vte(i, 2*mp1-1, k) = vte(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                    wto(i, 2*mp1-2, k) = wto(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                    wto(i, 2*mp1-1, k) = wto(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                end do
                                if (mlat /= 0) then
                                    vte(imid, 2*mp1-2, k) = vte(imid, 2*mp1-2, k) &
                                        -ci(mp1, np1, k)*wb(imid, mn)
                                    vte(imid, 2*mp1-1, k) = vte(imid, 2*mp1-1, k) &
                                        +cr(mp1, np1, k)*wb(imid, mn)
                                end if
                            end do
                        end do
                    end if
                end do
        end select vector_symmetry_cases

        do k=1, nt
            call util%hfft%backward(idv, nlon, vte(:, :, k), idv, wrfft)
            call util%hfft%backward(idv, nlon, wte(:, :, k), idv, wrfft)
        end do

        select case (ityp)
            case(0:1)
                do k=1, nt
                    do j=1, nlon
                        do i=1, imm1
                            vt(i, j, k) = HALF *(vte(i, j, k)+vto(i, j, k))
                            wt(i, j, k) = HALF *(wte(i, j, k)+wto(i, j, k))
                            vt(nlp1-i, j, k) = HALF *(vte(i, j, k)-vto(i, j, k))
                            wt(nlp1-i, j, k) = HALF *(wte(i, j, k)-wto(i, j, k))
                        end do
                    end do
                end do
            case default
                do k=1, nt
                    do j=1, nlon
                        do i=1, imm1
                            vt(i, j, k) = HALF *vte(i, j, k)
                            wt(i, j, k) = HALF *wte(i, j, k)
                        end do
                    end do
                end do
        end select

        if (mlat /= 0) then
            do k=1, nt
                do j=1, nlon
                    vt(imid, j, k) = HALF *vte(imid, j, k)
                    wt(imid, j, k) = HALF *wte(imid, j, k)
                end do
            end do
        end if

    end subroutine vtses_lower_utility_routine

    subroutine vtsesi_lower_utility_routine(nlat, nlon, imid, vb, wb, idz, vin, wzvin, dwork)
        integer(ip) :: i3
        integer(ip) :: idz
        integer(ip) :: imid
        integer(ip) :: m, mmax
        integer(ip) :: mn
        integer(ip) :: mp1
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: np1
        real(wp) :: vb(imid, *), wb(imid, *), vin(imid, nlat, 3), wzvin(*)
        real(wp) :: dwork(*)

        type(SpherepackUtility) :: util

        mmax = min(nlat, (nlon + 1)/2)

        call util%initialize_polar_components_regular_colat_deriv(nlat, nlon, wzvin, dwork)

        do mp1=1, mmax
            m = mp1-1
            call util%compute_polar_component(0, nlat, nlon, m, vin, i3, wzvin)
            do np1=mp1, nlat
                mn = m*(nlat-1)-(m*(m-1))/2+np1
                vb(1:imid, mn) = vin(1:imid, np1, i3)
            end do
        end do

        call util%initialize_azimuthal_components_regular_colat_deriv(nlat, nlon, wzvin, dwork)

        do mp1=1, mmax
            m = mp1-1
            call util%compute_azimuthal_component(0, nlat, nlon, m, vin, i3, wzvin)
            do np1=mp1, nlat
                mn = m*(nlat-1)-(m*(m-1))/2+np1
                wb(1:imid, mn) = vin(1:imid, np1, i3)
            end do
        end do

    end subroutine vtsesi_lower_utility_routine

end submodule colatitudinal_derivative_regular_grid_saved
