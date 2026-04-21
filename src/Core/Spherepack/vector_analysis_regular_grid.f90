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
submodule(vector_analysis_routines) vector_analysis_regular_grid

contains

    ! Purpose:
    !
    !     subroutine vhaec(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
    !                      mdab, ndab, wvhaec, ierror)
    !
    !     subroutine vhaec performs the vector spherical harmonic analysis
    !     on the vector field (v, w) and stores the result in the arrays
    !     br, bi, cr, and ci. v(i, j) and w(i, j) are the colatitudinal
    !     (measured from the north pole) and east longitudinal components
    !     respectively, located at colatitude theta(i) = (i-1)*pi/(nlat-1)
    !     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral
    !     representation of (v, w) is given at output parameters v, w in
    !     subroutine vhsec.
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
    !            nlon = 72 for a five degree grid. nlon must be greater
    !            than zero. the axisymmetric case corresponds to nlon=1.
    !            the efficiency of the computation is improved when nlon
    !            is a product of small prime numbers.
    !
    !     ityp   = 0  no symmetries exist about the equator. the analysis
    !                 is performed on the entire sphere.  i.e. on the
    !                 arrays v(i, j), w(i, j) for i=1, ..., nlat and
    !                 j=1, ..., nlon.
    !
    !            = 1  no symmetries exist about the equator. the analysis
    !                 is performed on the entire sphere.  i.e. on the
    !                 arrays v(i, j), w(i, j) for i=1, ..., nlat and
    !                 j=1, ..., nlon. the curl of (v, w) is zero. that is, 
    !                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0.
    !                 the coefficients cr and ci are zero.
    !
    !            = 2  no symmetries exist about the equator. the analysis
    !                 is performed on the entire sphere.  i.e. on the
    !                 arrays v(i, j), w(i, j) for i=1, ..., nlat and
    !                 j=1, ..., nlon. the divergence of (v, w) is zero. i.e., 
    !                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0.
    !                 the coefficients br and bi are zero.
    !
    !            = 3  v is symmetric and w is antisymmetric about the
    !                 equator. the analysis is performed on the northern
    !                 hemisphere only.  i.e., if nlat is odd the analysis
    !                 is performed on the arrays v(i, j), w(i, j) for
    !                 i=1, ..., (nlat + 1)/2 and j=1, ..., nlon. if nlat is
    !                 even the analysis is performed on the the arrays
    !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !            = 4  v is symmetric and w is antisymmetric about the
    !                 equator. the analysis is performed on the northern
    !                 hemisphere only.  i.e., if nlat is odd the analysis
    !                 is performed on the arrays v(i, j), w(i, j) for
    !                 i=1, ..., (nlat + 1)/2 and j=1, ..., nlon. if nlat is
    !                 even the analysis is performed on the the arrays
    !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !                 the curl of (v, w) is zero. that is, 
    !                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0.
    !                 the coefficients cr and ci are zero.
    !
    !            = 5  v is symmetric and w is antisymmetric about the
    !                 equator. the analysis is performed on the northern
    !                 hemisphere only.  i.e., if nlat is odd the analysis
    !                 is performed on the arrays v(i, j), w(i, j) for
    !                 i=1, ..., (nlat + 1)/2 and j=1, ..., nlon. if nlat is
    !                 even the analysis is performed on the the arrays
    !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !                 the divergence of (v, w) is zero. i.e., 
    !                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0.
    !                 the coefficients br and bi are zero.
    !
    !            = 6  v is antisymmetric and w is symmetric about the
    !                 equator. the analysis is performed on the northern
    !                 hemisphere only.  i.e., if nlat is odd the analysis
    !                 is performed on the arrays v(i, j), w(i, j) for
    !                 i=1, ..., (nlat + 1)/2 and j=1, ..., nlon. if nlat is
    !                 even the analysis is performed on the the arrays
    !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !            = 7  v is antisymmetric and w is symmetric about the
    !                 equator. the analysis is performed on the northern
    !                 hemisphere only.  i.e., if nlat is odd the analysis
    !                 is performed on the arrays v(i, j), w(i, j) for
    !                 i=1, ..., (nlat + 1)/2 and j=1, ..., nlon. if nlat is
    !                 even the analysis is performed on the the arrays
    !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !                 the curl of (v, w) is zero. that is, 
    !                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0.
    !                 the coefficients cr and ci are zero.
    !
    !            = 8  v is antisymmetric and w is symmetric about the
    !                 equator. the analysis is performed on the northern
    !                 hemisphere only.  i.e., if nlat is odd the analysis
    !                 is performed on the arrays v(i, j), w(i, j) for
    !                 i=1, ..., (nlat + 1)/2 and j=1, ..., nlon. if nlat is
    !                 even the analysis is performed on the the arrays
    !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !                 the divergence of (v, w) is zero. i.e., 
    !                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0.
    !                 the coefficients br and bi are zero.
    !
    !
    !     nt     the number of analyses.  in the program that calls vhaec, 
    !            the arrays v, w, br, bi, cr, and ci can be three dimensional
    !            in which case multiple analyses will be performed.
    !            the third index is the analysis index which assumes the
    !            values k=1, ..., nt.  for a single analysis set nt=1. the
    !            discription of the remaining parameters is simplified
    !            by assuming that nt=1 or that all the arrays are two
    !            dimensional.
    !
    !     v, w    two or three dimensional arrays (see input parameter nt)
    !            that contain the vector function to be analyzed.
    !            v is the colatitudnal component and w is the east
    !            longitudinal component. v(i, j), w(i, j) contain the
    !            components at colatitude theta(i) = (i-1)*pi/(nlat-1)
    !            and longitude phi(j) = (j-1)*2*pi/nlon. the index ranges
    !            are defined above at the input parameter ityp.
    !
    !     idvw   the first dimension of the arrays v, w as it appears in
    !            the program that calls vhaec. if ityp <= 2 then idvw
    !            must be at least nlat.  if ityp > 2 and nlat is
    !            even then idvw must be at least nlat/2. if ityp > 2
    !            and nlat is odd then idvw must be at least (nlat + 1)/2.
    !
    !     jdvw   the second dimension of the arrays v, w as it appears in
    !            the program that calls vhaec. jdvw must be at least nlon.
    !
    !     mdab   the first dimension of the arrays br, bi, cr, and ci as it
    !            appears in the program that calls vhaec. mdab must be at
    !            least min(nlat, nlon/2) if nlon is even or at least
    !            min(nlat, (nlon + 1)/2) if nlon is odd.
    !
    !     ndab   the second dimension of the arrays br, bi, cr, and ci as it
    !            appears in the program that calls vhaec. ndab must be at
    !            least nlat.
    !
    !     wvhaec an array which must be initialized by subroutine vhaeci.
    !            once initialized, wvhaec can be used repeatedly by vhaec
    !            as long as nlon and nlat remain unchanged.  wvhaec must
    !            not be altered between calls of vhaec.
    !
    !     lvhaec the dimension of the array wvhaec as it appears in the
    !            program that calls vhaec. define
    !
    !               l1 = min(nlat, nlon/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lvhaec must be at least
    !
    !            4*nlat*l2+3*max(l1-2, 0)*(2*nlat-l1-1)+nlon+15
    !
    !     output parameters
    !
    !     br, bi  two or three dimensional arrays (see input parameter nt)
    !     cr, ci  that contain the vector spherical harmonic coefficients
    !            in the spectral representation of v(i, j) and w(i, j) given
    !            in the discription of subroutine vhsec. br(mp1, np1), 
    !            bi(mp1, np1), cr(mp1, np1), and ci(mp1, np1) are computed
    !            for mp1=1, ..., mmax and np1=mp1, ..., nlat except for np1=nlat
    !            and odd mp1. mmax=min(nlat, nlon/2) if nlon is even or
    !            mmax=min(nlat, (nlon + 1)/2) if nlon is odd.
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
    !            = 9  error in the specification of lvhaec
    !
    module subroutine vhaec(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
        mdab, ndab, wvhaec, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: ityp
        integer(ip), intent(in)  :: nt
        real(wp),    intent(in)  :: v(idvw, jdvw, nt)
        real(wp),    intent(in)  :: w(idvw, jdvw, nt)
        integer(ip), intent(in)  :: idvw
        integer(ip), intent(in)  :: jdvw
        real(wp),    intent(out) :: br(mdab, ndab, nt)
        real(wp),    intent(out) :: bi(mdab, ndab, nt)
        real(wp),    intent(out) :: cr(mdab, ndab, nt)
        real(wp),    intent(out) :: ci(mdab, ndab, nt)
        integer(ip), intent(in)  :: mdab
        integer(ip), intent(in)  :: ndab
        real(wp),    intent(in)  :: wvhaec(:)
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip) :: mmax, lzz1
        integer(ip) :: labc, lnl, imid
        integer(ip) :: ist, idv, lwork
        integer(ip) :: required_wavetable_size
        type(VectorAnalysisUtility) :: util

        imid = (nlat + 1)/2
        mmax = min(nlat, (nlon + 1)/2)
        lzz1 = 2*nlat*imid
        labc = 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
        required_wavetable_size = 2 * (lzz1 + labc) + nlon + 15

         ! Check calling arguments
        call util%check_vector_analysis_inputs(nlat, nlon, ityp, idvw, jdvw, &
            mdab, ndab, nt, required_wavetable_size, wvhaec, ierror)

        ! Check error flag
        if (ierror /= 0) return

        ! Set required workspace size
        select case (ityp)
            case (0:2)
                lwork = nlat*(2*nt*nlon+max(6*imid, nlon))
            case default
                lwork = imid*(2*nt*nlon+max(6*nlat, nlon))
        end select

        select case (ityp)
            case (0:2)
                idv = nlat
                ist = imid
            case default
                idv = imid
                ist = 0
        end select

        lnl = nt*idv*nlon

        block
            real (wp)   :: work(lwork)
            integer(ip) :: iw1, iw2, iw3, iw4, iw5
            integer(ip) :: lwzvin, jw1, jw2

            ! Set workspace index pointers
            iw1 = ist+1
            iw2 = lnl+1
            iw3 = iw2+ist
            iw4 = iw2+lnl
            iw5 = iw4+3*imid*nlat

            ! Set wavetable index pointers
            lwzvin = lzz1+labc
            jw1 = lwzvin+1
            jw2 = jw1+lwzvin

            call vhaec_lower_utility_routine(nlat, nlon, ityp, nt, imid, idvw, jdvw, &
                v, w, mdab, ndab, br, bi, cr, ci, idv, work, work(iw1:), work(iw2:), work(iw3:), &
                work(iw4:), work(iw5:), wvhaec, wvhaec(jw1:), wvhaec(jw2:))
        end block

    end subroutine vhaec

    ! Purpose:
    !
    !     subroutine vhaeci(nlat, nlon, wvhaec, lvhaec, dwork, ldwork, ierror)
    !
    !     subroutine vhaeci initializes the array wvhaec which can then be
    !     used repeatedly by subroutine vhaec until nlat or nlon is changed.
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
    !            nlon = 72 for a five degree grid. nlon must be greater
    !            than zero. the axisymmetric case corresponds to nlon=1.
    !            the efficiency of the computation is improved when nlon
    !            is a product of small prime numbers.
    !
    !     lvhaec the dimension of the array wvhaec as it appears in the
    !            program that calls vhaec. define
    !
    !               l1 = min(nlat, nlon/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lvhaec must be at least
    !
    !            4*nlat*l2+3*max(l1-2, 0)*(2*nlat-l1-1)+nlon+15
    !
    !     output parameters
    !
    !     wvhaec an array which is initialized for use by subroutine vhaec.
    !            once initialized, wvhaec can be used repeatedly by vhaec
    !            as long as nlat or nlon remain unchanged.  wvhaec must not
    !            be altered between calls of vhaec.
    !
    !
    !     ierror = 0  no errors
    !            = 1  error in the specification of nlat
    !            = 2  error in the specification of nlon
    !            = 3  error in the specification of lvhaec
    !
    module subroutine vhaeci(nlat, nlon, wvhaec, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(out) :: wvhaec(:)
        integer(ip), intent(out) :: ierror

        ! Local variables
        type(SpherepackUtility) :: util
        integer(ip) :: imid
        integer(ip) :: iw1
        integer(ip) :: iw2
        integer(ip) :: labc
        integer(ip) :: lwzvin
        integer(ip) :: lzz1
        integer(ip) :: mmax, ldwork

        associate (lvhaec => size(wvhaec))

            imid = (nlat + 1)/2
            lzz1 = 2*nlat*imid
            mmax = min(nlat, (nlon + 1)/2)
            labc = 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2

            ! Check error flag
            if (nlat < 3) then
                ierror = 1
            else if (nlon < 1) then
                ierror = 2
            else if (lvhaec < 2*(lzz1+labc)+nlon+15) then
                ierror = 3
            else
                ierror = 0
            end if

            ! Check error flag
            if (ierror /= 0) return

            ! Set required workspace size
            ldwork = 2*nlat+2

            block
                real(wp) :: dwork(ldwork)

                call util%zvinit(nlat, nlon, wvhaec, dwork)

                ! Set workspace index pointers
                lwzvin = lzz1+labc
                iw1 = lwzvin+1
                iw2 = iw1+lwzvin

                call util%zwinit(nlat, nlon, wvhaec(iw1:), dwork)

                call util%hfft%initialize(nlon, wvhaec(iw2:))
            end block
        end associate

    end subroutine vhaeci

    subroutine vhaec_lower_utility_routine(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
        ndab, br, bi, cr, ci, idv, ve, vo, we, wo, zv, zw, wzvin, wzwin, wrfft)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: ityp
        integer(ip), intent(in)  :: nt
        integer(ip), intent(in)  :: imid
        integer(ip), intent(in)  :: idvw
        integer(ip), intent(in)  :: jdvw
        real(wp),    intent(in)  :: v(idvw, jdvw, nt)
        real(wp),    intent(in)  :: w(idvw, jdvw, nt)
        integer(ip), intent(in)  :: mdab
        integer(ip), intent(in)  :: ndab
        real(wp),    intent(out) :: br(mdab, ndab, nt)
        real(wp),    intent(out) :: bi(mdab, ndab, nt)
        real(wp),    intent(out) :: cr(mdab, ndab, nt)
        real(wp),    intent(out) :: ci(mdab, ndab, nt)
        integer(ip), intent(in)  :: idv
        real(wp),    intent(out) :: ve(idv, nlon, nt)
        real(wp),    intent(out) :: vo(idv, nlon, nt)
        real(wp),    intent(out) :: we(idv, nlon, nt)
        real(wp),    intent(out) :: wo(idv, nlon, nt)
        real(wp),    intent(out) :: zv(imid, nlat, 3)
        real(wp),    intent(out) :: zw(imid, nlat, 3)
        real(wp),    intent(in)  :: wzvin(:)
        real(wp),    intent(in)  :: wzwin(:)
        real(wp),    intent(in)  :: wrfft(:)

        ! Local variables
        integer(ip) :: i, imm1, iv, iw, k
        integer(ip) :: m, mmax, mp1, mp2
        integer(ip) :: odd_stride
        integer(ip) :: even_stride, np1
        type(VectorAnalysisUtility) :: util

        call util%analysis_setup(idvw, jdvw, mdab, ndab, &
            bi, br, ci, cr, even_stride, idv, imid, imm1, ityp, &
            mmax, nlat, nlon, nt, odd_stride,  v, ve, vo, w, we, wo, wrfft)

        select case (ityp)
            case(0)

                ! case ityp=0 ,  no symmetries
                call util%zvin(0, nlat, nlon, 0, zv, iv, wzvin)

                ! case m=0
                do k=1, nt
                    do i=1, imid
                        do np1=2, even_stride, 2
                            br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*ve(i, 1, k)
                            cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*we(i, 1, k)
                        end do
                    end do
                end do
                do k=1, nt
                    do i=1, imm1
                        do np1=3, odd_stride, 2
                            br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*vo(i, 1, k)
                            cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*wo(i, 1, k)
                        end do
                    end do
                end do

                ! case m = 1 through nlat-1
                if (mmax < 2) return

                main_loop_case_0: do mp1=2, mmax
                    m = mp1-1
                    mp2 = mp1+1

                    call util%zvin(0, nlat, nlon, m, zv, iv, wzvin)
                    call util%zwin(0, nlat, nlon, m, zw, iw, wzwin)

                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do i=1, imm1
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*vo(i, 2*mp1-2, k) &
                                        +zw(i, np1, iw)*we(i, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*vo(i, 2*mp1-1, k) &
                                        -zw(i, np1, iw)*we(i, 2*mp1-2, k)
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-2, k) &
                                        +zw(i, np1, iw)*ve(i, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-1, k) &
                                        -zw(i, np1, iw)*ve(i, 2*mp1-2, k)
                                end do
                            end do
                        end do
                        if (odd(nlat)) then
                            do k=1, nt
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+zw(imid, np1, iw)*we(imid, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)-zw(imid, np1, iw)*we(imid, 2*mp1-2, k)
                                    cr(mp1, np1, k) = cr(mp1, np1, k)+zw(imid, np1, iw)*ve(imid, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-zw(imid, np1, iw)*ve(imid, 2*mp1-2, k)
                                end do
                            end do
                        end if
                    end if

                    if (mp2 > even_stride) cycle main_loop_case_0

                    do k=1, nt
                        do i=1, imm1
                            do np1=mp2, even_stride, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-2, k) &
                                    +zw(i, np1, iw)*wo(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-1, k) &
                                    -zw(i, np1, iw)*wo(i, 2*mp1-2, k)
                                cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*we(i, 2*mp1-2, k) &
                                    +zw(i, np1, iw)*vo(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*we(i, 2*mp1-1, k) &
                                    -zw(i, np1, iw)*vo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do

                    if (even(nlat)) cycle main_loop_case_0

                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-2, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-1, k)
                            cr(mp1, np1, k) = cr(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-2, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-1, k)
                        end do
                    end do
                end do main_loop_case_0
            case (1)
                                !
                                ! case ityp=1 ,  no symmetries but cr and ci equal zero
                                !
                call util%zvin(0, nlat, nlon, 0, zv, iv, wzvin)
                !
                ! case m=0
                !
                do k=1, nt
                    do i=1, imid
                        do np1=2, even_stride, 2
                            br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*ve(i, 1, k)
                        end do
                    end do
                end do

                do k=1, nt
                    do i=1, imm1
                        do np1=3, odd_stride, 2
                            br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*vo(i, 1, k)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) return
                main_loop_case_1: do mp1=2, mmax
                    m = mp1-1
                    mp2 = mp1+1
                    call util%zvin(0, nlat, nlon, m, zv, iv, wzvin)
                    call util%zwin(0, nlat, nlon, m, zw, iw, wzwin)
                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do i=1, imm1
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*vo(i, 2*mp1-2, k) &
                                        +zw(i, np1, iw)*we(i, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*vo(i, 2*mp1-1, k) &
                                        -zw(i, np1, iw)*we(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (odd(nlat)) then
                            do k=1, nt
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+zw(imid, np1, iw)*we(imid, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)-zw(imid, np1, iw)*we(imid, 2*mp1-2, k)
                                end do
                            end do
                        end if
                    end if
                    if (mp2 > even_stride) cycle main_loop_case_1
                    do k=1, nt
                        do i=1, imm1
                            do np1=mp2, even_stride, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-2, k) &
                                    +zw(i, np1, iw)*wo(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-1, k) &
                                    -zw(i, np1, iw)*wo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do

                    if (even(nlat)) cycle main_loop_case_1
                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-2, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-1, k)
                        end do
                    end do
                end do main_loop_case_1
            case (2)
                               !
                               ! case ityp=2 ,  no symmetries but br and bi equal zero
                               !
                call util%zvin(0, nlat, nlon, 0, zv, iv, wzvin)
                !
                ! case m=0
                !
                do k=1, nt
                    do i=1, imid
                        do np1=2, even_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*we(i, 1, k)
                        end do
                    end do
                end do

                do k=1, nt
                    do i=1, imm1
                        do np1=3, odd_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*wo(i, 1, k)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) return
                main_loop_case_2: do mp1=2, mmax
                    m = mp1-1
                    mp2 = mp1+1
                    call util%zvin(0, nlat, nlon, m, zv, iv, wzvin)
                    call util%zwin(0, nlat, nlon, m, zw, iw, wzwin)
                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do i=1, imm1
                                do np1=mp1, odd_stride, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-2, k) &
                                        +zw(i, np1, iw)*ve(i, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-1, k) &
                                        -zw(i, np1, iw)*ve(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (odd(nlat)) then
                            do k=1, nt
                                do np1=mp1, odd_stride, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)+zw(imid, np1, iw)*ve(imid, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-zw(imid, np1, iw)*ve(imid, 2*mp1-2, k)
                                end do
                            end do
                        end if
                    end if

                    if (mp2 > even_stride) cycle main_loop_case_2
                    do k=1, nt
                        do i=1, imm1
                            do np1=mp2, even_stride, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*we(i, 2*mp1-2, k) &
                                    +zw(i, np1, iw)*vo(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*we(i, 2*mp1-1, k) &
                                    -zw(i, np1, iw)*vo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do
                    if (even(nlat)) cycle main_loop_case_2
                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-2, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-1, k)
                        end do
                    end do
                end do main_loop_case_2
            case (3)
                                !
                                ! case ityp=3 ,  v even , w odd
                                !
                call util%zvin(0, nlat, nlon, 0, zv, iv, wzvin)
                !
                ! case m=0
                !
                do k=1, nt
                    do i=1, imid
                        do np1=2, even_stride, 2
                            br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*ve(i, 1, k)
                        end do
                    end do
                end do
                do k=1, nt
                    do i=1, imm1
                        do np1=3, odd_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*wo(i, 1, k)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) return
                main_loop_case_3: do mp1=2, mmax
                    m = mp1-1
                    mp2 = mp1+1
                    call util%zvin(0, nlat, nlon, m, zv, iv, wzvin)
                    call util%zwin(0, nlat, nlon, m, zw, iw, wzwin)
                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do i=1, imm1
                                do np1=mp1, odd_stride, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-2, k) &
                                        +zw(i, np1, iw)*ve(i, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-1, k) &
                                        -zw(i, np1, iw)*ve(i, 2*mp1-2, k)
                                end do
                            end do
                        end do
                        if (odd(nlat)) then
                            do k=1, nt
                                do np1=mp1, odd_stride, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)+zw(imid, np1, iw)*ve(imid, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-zw(imid, np1, iw)*ve(imid, 2*mp1-2, k)
                                end do
                            end do
                        end if
                    end if
                    if (mp2 > even_stride) cycle main_loop_case_3
                    do k=1, nt
                        do i=1, imm1
                            do np1=mp2, even_stride, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-2, k) &
                                    +zw(i, np1, iw)*wo(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-1, k) &
                                    -zw(i, np1, iw)*wo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do
                    if (even(nlat)) cycle main_loop_case_3
                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-2, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-1, k)
                        end do
                    end do
                end do main_loop_case_3
            case (4)
                                !
                                ! case ityp=4 ,  v even, w odd, and cr and ci equal 0.
                                !
                call util%zvin(1, nlat, nlon, 0, zv, iv, wzvin)
                !
                ! case m=0
                !
                do k=1, nt
                    do i=1, imid
                        do np1=2, even_stride, 2
                            br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*ve(i, 1, k)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) return
                main_loop_case_4: do mp1=2, mmax
                    m = mp1-1
                    mp2 = mp1+1
                    call util%zvin(1, nlat, nlon, m, zv, iv, wzvin)
                    call util%zwin(1, nlat, nlon, m, zw, iw, wzwin)
                    if (mp2 > even_stride) cycle main_loop_case_4
                    do k=1, nt
                        do i=1, imm1
                            do np1=mp2, even_stride, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-2, k) &
                                    +zw(i, np1, iw)*wo(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-1, k) &
                                    -zw(i, np1, iw)*wo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do
                    if (even(nlat)) cycle main_loop_case_4
                    do k=1, nt
                        do  np1=mp2, even_stride, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-2, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-1, k)
                        end do
                    end do
                end do main_loop_case_4
            case (5)
                                !
                                ! case ityp=5   v even, w odd, and br and bi equal zero
                                !
                call util%zvin(2, nlat, nlon, 0, zv, iv, wzvin)
                !
                ! case m=0
                !
                do k=1, nt
                    do i=1, imm1
                        do np1=3, odd_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*wo(i, 1, k)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) return
                main_loop_case_5: do mp1=2, mmax
                    m = mp1-1
                    mp2 = mp1+1
                    call util%zvin(2, nlat, nlon, m, zv, iv, wzvin)
                    call util%zwin(2, nlat, nlon, m, zw, iw, wzwin)
                    if (mp1 > odd_stride) cycle main_loop_case_5
                    do k=1, nt
                        do i=1, imm1
                            do np1=mp1, odd_stride, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-2, k) &
                                    +zw(i, np1, iw)*ve(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-1, k) &
                                    -zw(i, np1, iw)*ve(i, 2*mp1-2, k)
                            end do
                        end do
                    end do
                    if (even(nlat)) cycle main_loop_case_5
                    do k=1, nt
                        do np1=mp1, odd_stride, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)+zw(imid, np1, iw)*ve(imid, 2*mp1-1, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-zw(imid, np1, iw)*ve(imid, 2*mp1-2, k)
                        end do
                    end do
                end do main_loop_case_5
            case (6)
                                !
                                ! case ityp=6 ,  v odd , w even
                                !
                call util%zvin(0, nlat, nlon, 0, zv, iv, wzvin)
                !
                ! case m=0
                !
                do k=1, nt
                    do i=1, imid
                        do np1=2, even_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*we(i, 1, k)
                        end do
                    end do
                end do
                do k=1, nt
                    do i=1, imm1
                        do np1=3, odd_stride, 2
                            br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*vo(i, 1, k)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) return
                main_loop_case_6: do mp1=2, mmax
                    m = mp1-1
                    mp2 = mp1+1
                    call util%zvin(0, nlat, nlon, m, zv, iv, wzvin)
                    call util%zwin(0, nlat, nlon, m, zw, iw, wzwin)
                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do i=1, imm1
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*vo(i, 2*mp1-2, k) &
                                        +zw(i, np1, iw)*we(i, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*vo(i, 2*mp1-1, k) &
                                        -zw(i, np1, iw)*we(i, 2*mp1-2, k)
                                end do
                            end do
                        end do
                        if (odd(nlat)) then
                            do k=1, nt
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+zw(imid, np1, iw)*we(imid, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)-zw(imid, np1, iw)*we(imid, 2*mp1-2, k)
                                end do
                            end do
                        end if
                    end if

                    if (mp2 > even_stride) cycle main_loop_case_6
                    do k=1, nt
                        do i=1, imm1
                            do np1=mp2, even_stride, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*we(i, 2*mp1-2, k) &
                                    +zw(i, np1, iw)*vo(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*we(i, 2*mp1-1, k) &
                                    -zw(i, np1, iw)*vo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do
                    if (even(nlat)) cycle main_loop_case_6
                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-2, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-1, k)
                        end do
                    end do
                end do main_loop_case_6
            case (7)

                ! case ityp=7   v odd, w even, and cr and ci equal zero
                call util%zvin(2, nlat, nlon, 0, zv, iv, wzvin)

                ! case m=0
                do k=1, nt
                    do i=1, imm1
                        do np1=3, odd_stride, 2
                            br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*vo(i, 1, k)
                        end do
                    end do
                end do

                ! case m = 1 through nlat-1
                if (mmax < 2) return
                main_loop_case_7: do mp1=2, mmax
                    m = mp1-1
                    mp2 = mp1+1
                    call util%zvin(2, nlat, nlon, m, zv, iv, wzvin)
                    call util%zwin(2, nlat, nlon, m, zw, iw, wzwin)
                    if (mp1 > odd_stride) cycle main_loop_case_7
                    do k=1, nt
                        do i=1, imm1
                            do np1=mp1, odd_stride, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*vo(i, 2*mp1-2, k) &
                                    +zw(i, np1, iw)*we(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*vo(i, 2*mp1-1, k) &
                                    -zw(i, np1, iw)*we(i, 2*mp1-2, k)
                            end do
                        end do
                    end do
                    if (even(nlat)) cycle main_loop_case_7
                    do k=1, nt
                        do np1=mp1, odd_stride, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+zw(imid, np1, iw)*we(imid, 2*mp1-1, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)-zw(imid, np1, iw)*we(imid, 2*mp1-2, k)
                        end do
                    end do
                end do main_loop_case_7
            case (8)

                ! case ityp=8   v odd, w even, and both br and bi equal zero
                call util%zvin(1, nlat, nlon, 0, zv, iv, wzvin)

                ! case m=0
                do k=1, nt
                    do i=1, imid
                        do np1=2, even_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*we(i, 1, k)
                        end do
                    end do
                end do

                ! case m = 1 through nlat-1
                if (mmax < 2) return

                main_loop_case_8: do mp1=2, mmax
                    m = mp1-1
                    mp2 = mp1+1
                    call util%zvin(1, nlat, nlon, m, zv, iv, wzvin)
                    call util%zwin(1, nlat, nlon, m, zw, iw, wzwin)
                    if (mp2 > even_stride) cycle main_loop_case_8
                    do k=1, nt
                        do i=1, imm1
                            do np1=mp2, even_stride, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*we(i, 2*mp1-2, k) &
                                    +zw(i, np1, iw)*vo(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*we(i, 2*mp1-1, k) &
                                    -zw(i, np1, iw)*vo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do
                    if (even(nlat)) cycle main_loop_case_8
                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-2, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-1, k)
                        end do
                    end do
                end do main_loop_case_8
        end select

    end subroutine vhaec_lower_utility_routine

end submodule vector_analysis_regular_grid
