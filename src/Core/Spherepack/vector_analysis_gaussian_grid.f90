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
submodule(vector_analysis_routines) vector_analysis_gaussian_grid

contains

    !     subroutine vhagc(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
    !                      mdab, ndab, wvhagc, ierror)
    !
    !     subroutine vhagc performs the vector spherical harmonic analysis
    !     on the vector field (v, w) and stores the result in the arrays
    !     br, bi, cr, and ci. v(i, j) and w(i, j) are the colatitudinal
    !     (measured from the north pole) and east longitudinal components
    !     respectively, located at the gaussian colatitude point theta(i)
    !     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral
    !     representation of (v, w) is given at output parameters v, w in
    !     subroutine vhsec.
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
    !     nt     the number of analyses.  in the program that calls vhagc,
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
    !            the program that calls vhagc. if ityp <= 2 then idvw
    !            must be at least nlat.  if ityp > 2 and nlat is
    !            even then idvw must be at least nlat/2. if ityp > 2
    !            and nlat is odd then idvw must be at least (nlat + 1)/2.
    !
    !     jdvw   the second dimension of the arrays v, w as it appears in
    !            the program that calls vhagc. jdvw must be at least nlon.
    !
    !     mdab   the first dimension of the arrays br, bi, cr, and ci as it
    !            appears in the program that calls vhagc. mdab must be at
    !            least min(nlat, nlon/2) if nlon is even or at least
    !            min(nlat, (nlon + 1)/2) if nlon is odd.
    !
    !     ndab   the second dimension of the arrays br, bi, cr, and ci as it
    !            appears in the program that calls vhagc. ndab must be at
    !            least nlat.
    !
    !     wvhagc an array which must be initialized by subroutine vhagci.
    !            once initialized, wvhagc can be used repeatedly by vhagc
    !            as long as nlon and nlat remain unchanged.  wvhagc must
    !            not be altered between calls of vhagc.
    !
    !     lvhagc the dimension of the array wvhagc as it appears in the
    !            program that calls vhagc. define
    !
    !               l1 = min(nlat, nlon/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lvhagc must be at least
    !
    !               4*nlat*l2+3*max(l1-2, 0)*(2*nlat-l1-1)+nlon+l2+15
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
    !            = 9  error in the specification of lvhagc
    !
    module subroutine vhagc(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
        mdab, ndab, wvhagc, ierror)

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
        real(wp),    intent(in)  :: wvhagc(:)
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip) :: imid, labc, lzz1
        integer(ip) :: required_wavetable_size, mmax, lwork
        type(VectorAnalysisUtility) :: util

        imid = (nlat + 1)/2
        mmax = min(nlat, (nlon + 1)/2)
        lzz1 = 2 * nlat * imid
        labc = 3 * (max(mmax - 2, 0) * ((2 * nlat) - mmax - 1))/2
        required_wavetable_size= 2 * (lzz1 + labc) + nlon + imid + 15

        ! Check calling arguments
        call util%check_vector_analysis_inputs(nlat, nlon, ityp, idvw, jdvw, &
            mdab, ndab, nt, required_wavetable_size, wvhagc, ierror)

        ! Check error flag
        if (ierror /= 0) return

        ! Set required workspace size
        select case (ityp)
            case (0:2)
                lwork = nlat * (4 * nlon * nt + (6 * imid))
            case default
                lwork = imid * (4 * nlon * nt + (6 * nlat))
        end select

        block
            integer(ip) :: idv, ist, lnl
            integer(ip) :: iw1, iw2, iw3, iw4, iw5
            integer(ip) :: jw1, jw2, jw3, lwzvin
            real(wp)    :: work(lwork)

            select case (ityp)
                case (0:2)
                    idv = nlat
                    ist = imid
                case default
                    idv = imid
                    ist = 0
            end select

            lnl = nt*idv*nlon

            ! Set workspace index pointers
            iw1 = ist+1
            iw2 = lnl+1
            iw3 = iw2+ist
            iw4 = iw2+lnl
            iw5 = iw4+3*imid*nlat

            ! Set wavetable index pointers
            lwzvin = lzz1+labc
            jw1 = (nlat + 1)/2+1
            jw2 = jw1+lwzvin
            jw3 = jw2+lwzvin

            call vhagc_lower_utility_routine(nlat, nlon, ityp, nt, imid, idvw, jdvw, &
                v, w, mdab, ndab, br, bi, cr, ci, idv, work, work(iw1:), work(iw2:), &
                work(iw3:), work(iw4:), work(iw5:), wvhagc, wvhagc(jw1:), wvhagc(jw2:), wvhagc(jw3:))
        end block

    end subroutine vhagc

    !     subroutine vhagci(nlat, nlon, wvhagc, ierror)
    !
    !     subroutine vhagci initializes the array wvhagc which can then be
    !     used repeatedly by subroutine vhagc until nlat or nlon is changed.
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
    !     lvhagc the dimension of the array wvhagc as it appears in the
    !            program that calls vhagci.  define
    !
    !               l1 = min(nlat, nlon/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lvhagc must be at least
    !
    !               4*nlat*l2+3*max(l1-2, 0)*(2*nlat-l1-1)+nlon+l2+15
    !
    !     output parameters
    !
    !     wvhagc an array which is initialized for use by subroutine vhagc.
    !            once initialized, wvhagc can be used repeatedly by vhagc
    !            as long as nlat and nlon remain unchanged.  wvhagc must not
    !            be altered between calls of vhagc.
    !
    !
    !     ierror = 0  no errors
    !            = 1  error in the specification of nlat
    !            = 2  error in the specification of nlon
    !            = 3  error in the specification of lvhagc
    !
    module subroutine vhagci(nlat, nlon, wvhagc, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(out) :: wvhagc(:)
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip) :: imid, labc, lzz1
        integer(ip) :: mmax, ldwork, required_wavetable_size
        type(SpherepackUtility) :: util

        imid = (nlat + 1)/2
        lzz1 = 2*nlat*imid
        mmax = min(nlat, (nlon + 1)/2)
        labc = 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
        imid = (nlat + 1)/2
        required_wavetable_size = 2*(lzz1+labc)+nlon+imid+15

         ! Check calling arguments
        if (nlat < 3) then
            ierror = 1
        else if (nlon < 1) then
            ierror = 2
        else if (size(wvhagc) < required_wavetable_size) then
            ierror = 3
        else
            ierror = 0
        end if

        ! Check error flag
        if (ierror /= 0) return

        ! Set required workspace size
        ldwork = (2 * nlat * (nlat + 1)) + 1

        block
            real(wp)    :: dwork(ldwork)
            integer(ip) :: jw1, jw2, jw3, iw1, iw2, iw3
            integer(ip) :: iwrk, lwvbin

            ! Compute gaussian points in first nlat+1 words of dwork
            jw1 = 1
            jw2 = jw1 + nlat
            jw3 = jw2 + nlat

            call compute_gaussian_latitudes_and_weights( &
                nlat, dwork(jw1:), dwork(jw2:), ierror)

            ! Set first imid words of real weights in dwork
            ! in first imid words of wvhagc
            imid = (nlat + 1)/2
            wvhagc(:imid) = dwork(jw2:imid)

            ! first nlat+1 words of dwork contain  double theta
            iwrk = imid + 1
            iw1 = iwrk
            lwvbin = lzz1 + labc
            iw2 = iw1 + lwvbin
            iw3 = iw2 + lwvbin

            call util%initialize_polar_components_gaussian_grid( &
                nlat, nlon, dwork, wvhagc(iw1:), dwork(iwrk:))

            call util%initialize_azimuthal_components_gaussian_grid( &
                nlat, nlon, dwork, wvhagc(iw2:), dwork(iwrk:))

            call util%hfft%initialize(nlon, wvhagc(iw3:))
        end block

    end subroutine vhagci

    subroutine vhagc_lower_utility_routine(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
        ndab, br, bi, cr, ci, idv, ve, vo, we, wo, vb, wb, wts, wvbin, wwbin, wrfft)

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
        
        real(wp) :: tv
        real(wp) :: tve1
        real(wp) :: tve2
        real(wp) :: tvo1
        real(wp) :: tvo2
        real(wp) :: tw
        real(wp) :: twe1
        real(wp) :: twe2
        real(wp) :: two1
        real(wp) :: two2
        real(wp) :: v
        real(wp) :: vb
        real(wp) :: ve
        real(wp) :: vo
        real(wp) :: w
        real(wp) :: wb
        real(wp) :: we
        real(wp) :: wo
        real(wp) :: wrfft(:)
        real(wp) :: wts(nlat)
        real(wp) :: wvbin
        real(wp) :: wwbin
        dimension v(idvw, jdvw, nt), w(idvw, jdvw, nt), br(mdab, ndab, nt), &
            bi(mdab, ndab, nt), cr(mdab, ndab, nt), ci(mdab, ndab, nt), &
            ve(idv, nlon, nt), vo(idv, nlon, nt), we(idv, nlon, nt), &
            wo(idv, nlon, nt), wvbin(*), wwbin(*), &
            vb(imid, nlat, 3), wb(imid, nlat, 3)

        type(VectorAnalysisUtility) :: util

        call util%analysis_setup(idvw, jdvw, mdab, ndab, &
            bi, br, ci, cr, even_stride, idv, imid, imm1, ityp, &
            mmax, nlat, nlon, nt, odd_stride,  v, ve, vo, w, we, wo, wrfft)

        vector_symmetry_cases: select case (ityp)
            case (0)

                ! Case ityp=0 ,  no symmetries

                call util%compute_polar_component(0, nlat, nlon, 0, vb, iv, wvbin)

                ! Case m=0
                do k=1, nt
                    do i=1, imid
                        tv = ve(i, 1, k) * wts(i)
                        tw = we(i, 1, k) * wts(i)
                        do np1=2, even_stride, 2
                            br(1, np1, k) = br(1, np1, k)+vb(i, np1, iv)*tv
                            cr(1, np1, k) = cr(1, np1, k)-vb(i, np1, iv)*tw
                        end do
                    end do
                end do

                do k=1, nt
                    do i=1, imm1
                        tv = vo(i, 1, k)*wts(i)
                        tw = wo(i, 1, k)*wts(i)
                        do np1=3, odd_stride, 2
                            br(1, np1, k) = br(1, np1, k)+vb(i, np1, iv)*tv
                            cr(1, np1, k) = cr(1, np1, k)-vb(i, np1, iv)*tw
                        end do
                    end do
                end do

                ! Case m = 1 through nlat-1
                if (mmax < 2) return

                main_loop_case_0: do mp1=2, mmax
                    m = mp1-1
                    mp2 = mp1+1
                    call util%compute_polar_component(0, nlat, nlon, m, vb, iv, wvbin)
                    call util%compute_azimuthal_component(0, nlat, nlon, m, wb, iw, wwbin)

                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do i=1, imm1

                                ! Set temps to optimize quadrature
                                tvo1 = vo(i, 2*mp1-1, k)*wts(i)
                                tvo2 = vo(i, 2*mp1-2, k)*wts(i)
                                tve1 = ve(i, 2*mp1-1, k)*wts(i)
                                tve2 = ve(i, 2*mp1-2, k)*wts(i)
                                two1 = wo(i, 2*mp1-1, k)*wts(i)
                                two2 = wo(i, 2*mp1-2, k)*wts(i)
                                twe1 = we(i, 2*mp1-1, k)*wts(i)
                                twe2 = we(i, 2*mp1-2, k)*wts(i)

                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*tvo2 &
                                        +wb(i, np1, iw)*twe1
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*tvo1 &
                                        -wb(i, np1, iw)*twe2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*two2 &
                                        +wb(i, np1, iw)*tve1
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*two1 &
                                        -wb(i, np1, iw)*tve2
                                end do
                            end do
                        end do
                        if (odd(nlat))then
                            i = imid
                            do k=1, nt
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k)=br(mp1, np1, k)+wb(i, np1, iw)*we(i, 2*mp1-1, k)*wts(i)
                                    bi(mp1, np1, k)=bi(mp1, np1, k)-wb(i, np1, iw)*we(i, 2*mp1-2, k)*wts(i)
                                    cr(mp1, np1, k)=cr(mp1, np1, k)+wb(i, np1, iw)*ve(i, 2*mp1-1, k)*wts(i)
                                    ci(mp1, np1, k)=ci(mp1, np1, k)-wb(i, np1, iw)*ve(i, 2*mp1-2, k)*wts(i)
                                end do
                            end do
                        end if
                    end if

                    if (mp2 > even_stride) cycle main_loop_case_0

                    do k=1, nt
                        do i=1, imm1
                            tvo1 = vo(i, 2*mp1-1, k)*wts(i)
                            tvo2 = vo(i, 2*mp1-2, k)*wts(i)
                            tve1 = ve(i, 2*mp1-1, k)*wts(i)
                            tve2 = ve(i, 2*mp1-2, k)*wts(i)
                            two1 = wo(i, 2*mp1-1, k)*wts(i)
                            two2 = wo(i, 2*mp1-2, k)*wts(i)
                            twe1 = we(i, 2*mp1-1, k)*wts(i)
                            twe2 = we(i, 2*mp1-2, k)*wts(i)
                            do np1=mp2, even_stride, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*tve2 &
                                    +wb(i, np1, iw)*two1
                                bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*tve1 &
                                    -wb(i, np1, iw)*two2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*twe2 &
                                    +wb(i, np1, iw)*tvo1
                                ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*twe1 &
                                    -wb(i, np1, iw)*tvo2
                            end do
                        end do
                    end do

                    if (even(nlat)) cycle main_loop_case_0

                    i = imid

                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            br(mp1, np1, k)=br(mp1, np1, k)+vb(i, np1, iv)*ve(i, 2*mp1-2, k)*wts(i)
                            bi(mp1, np1, k)=bi(mp1, np1, k)+vb(i, np1, iv)*ve(i, 2*mp1-1, k)*wts(i)
                            cr(mp1, np1, k)=cr(mp1, np1, k)-vb(i, np1, iv)*we(i, 2*mp1-2, k)*wts(i)
                            ci(mp1, np1, k)=ci(mp1, np1, k)-vb(i, np1, iv)*we(i, 2*mp1-1, k)*wts(i)
                        end do
                    end do
                end do main_loop_case_0
            case(1)
                !
                ! case ityp=1 ,  no symmetries but cr and ci equal zero
                !
                call util%compute_polar_component(0, nlat, nlon, 0, vb, iv, wvbin)
                !
                ! case m=0
                !
                do k=1, nt
                    do i=1, imid
                        tv = ve(i, 1, k)*wts(i)
                        do np1=2, even_stride, 2
                            br(1, np1, k) = br(1, np1, k)+vb(i, np1, iv)*tv
                        end do
                    end do
                end do

                do k=1, nt
                    do i=1, imm1
                        tv = vo(i, 1, k)*wts(i)
                        do np1=3, odd_stride, 2
                            br(1, np1, k) = br(1, np1, k)+vb(i, np1, iv)*tv
                        end do
                    end do
                end do

                ! case m = 1 through nlat-1
                if (mmax < 2) return

                main_loop_case_1: do mp1=2, mmax

                    m = mp1-1
                    mp2 = mp1+1

                    call util%compute_polar_component(0, nlat, nlon, m, vb, iv, wvbin)
                    call util%compute_azimuthal_component(0, nlat, nlon, m, wb, iw, wwbin)

                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do i=1, imm1
                                tvo1 = vo(i, 2*mp1-1, k)*wts(i)
                                tvo2 = vo(i, 2*mp1-2, k)*wts(i)
                                twe1 = we(i, 2*mp1-1, k)*wts(i)
                                twe2 = we(i, 2*mp1-2, k)*wts(i)
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*tvo2 &
                                        +wb(i, np1, iw)*twe1
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*tvo1 &
                                        -wb(i, np1, iw)*twe2
                                end do
                            end do
                        end do
                        if (odd(nlat))then
                            i = imid
                            do k=1, nt
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+wb(i, np1, iw)*we(i, 2*mp1-1, k)*wts(i)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)-wb(i, np1, iw)*we(i, 2*mp1-2, k)*wts(i)
                                end do
                            end do
                        end if
                    end if

                    if (mp2 > even_stride) cycle main_loop_case_1

                    do k=1, nt
                        do i=1, imm1
                            tve1 = ve(i, 2*mp1-1, k)*wts(i)
                            tve2 = ve(i, 2*mp1-2, k)*wts(i)
                            two1 = wo(i, 2*mp1-1, k)*wts(i)
                            two2 = wo(i, 2*mp1-2, k)*wts(i)
                            do np1=mp2, even_stride, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*tve2 &
                                    +wb(i, np1, iw)*two1
                                bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*tve1 &
                                    -wb(i, np1, iw)*two2
                            end do
                        end do
                    end do

                    if (even(nlat)) cycle main_loop_case_1

                    i = imid

                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*ve(i, 2*mp1-2, k)*wts(i)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*ve(i, 2*mp1-1, k)*wts(i)
                        end do
                    end do
                end do main_loop_case_1
            case(2)
                !
                ! case ityp=2 ,  no symmetries but br and bi equal zero
                !
                call util%compute_polar_component(0, nlat, nlon, 0, vb, iv, wvbin)
                !
                ! case m=0
                !
                do k=1, nt
                    do i=1, imid
                        tw = we(i, 1, k)*wts(i)
                        do np1=2, even_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-vb(i, np1, iv)*tw
                        end do
                    end do
                end do

                do k=1, nt
                    do i=1, imm1
                        tw = wo(i, 1, k)*wts(i)
                        do np1=3, odd_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-vb(i, np1, iv)*tw
                        end do
                    end do
                end do

                ! case m = 1 through nlat-1
                if (mmax < 2) return

                main_loop_case_2: do mp1=2, mmax

                    m = mp1-1
                    mp2 = mp1+1

                    call util%compute_polar_component(0, nlat, nlon, m, vb, iv, wvbin)
                    call util%compute_azimuthal_component(0, nlat, nlon, m, wb, iw, wwbin)

                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do i=1, imm1
                                tve1 = ve(i, 2*mp1-1, k)*wts(i)
                                tve2 = ve(i, 2*mp1-2, k)*wts(i)
                                two1 = wo(i, 2*mp1-1, k)*wts(i)
                                two2 = wo(i, 2*mp1-2, k)*wts(i)
                                do np1=mp1, odd_stride, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*two2 &
                                        +wb(i, np1, iw)*tve1
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*two1 &
                                        -wb(i, np1, iw)*tve2
                                end do
                            end do
                        end do
                        if (odd(nlat))then
                            i = imid
                            do k=1, nt
                                do np1=mp1, odd_stride, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)+wb(i, np1, iw)*ve(i, 2*mp1-1, k)*wts(i)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-wb(i, np1, iw)*ve(i, 2*mp1-2, k)*wts(i)
                                end do
                            end do
                        end if
                    end if

                    if (mp2 > even_stride) cycle main_loop_case_2

                    do k=1, nt
                        do i=1, imm1
                            twe1 = we(i, 2*mp1-1, k)*wts(i)
                            twe2 = we(i, 2*mp1-2, k)*wts(i)
                            tvo1 = vo(i, 2*mp1-1, k)*wts(i)
                            tvo2 = vo(i, 2*mp1-2, k)*wts(i)
                            do np1=mp2, even_stride, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*twe2 &
                                    +wb(i, np1, iw)*tvo1
                                ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*twe1 &
                                    -wb(i, np1, iw)*tvo2
                            end do
                        end do
                    end do

                    if (even(nlat)) cycle main_loop_case_2

                    i = imid

                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*we(i, 2*mp1-2, k)*wts(i)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*we(i, 2*mp1-1, k)*wts(i)
                        end do
                    end do
                end do main_loop_case_2
            case(3)
                !
                ! case ityp=3 ,  v even , w odd
                !
                call util%compute_polar_component(0, nlat, nlon, 0, vb, iv, wvbin)
                !
                ! case m=0
                !
                do k=1, nt
                    do i=1, imid
                        tv = ve(i, 1, k)*wts(i)
                        do np1=2, even_stride, 2
                            br(1, np1, k) = br(1, np1, k)+vb(i, np1, iv)*tv
                        end do
                    end do
                end do

                do k=1, nt
                    do i=1, imm1
                        tw = wo(i, 1, k)*wts(i)
                        do np1=3, odd_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-vb(i, np1, iv)*tw
                        end do
                    end do
                end do

                ! case m = 1 through nlat-1
                if (mmax < 2) return

                main_loop_case_3: do mp1=2, mmax

                    m = mp1-1
                    mp2 = mp1+1

                    call util%compute_polar_component(0, nlat, nlon, m, vb, iv, wvbin)
                    call util%compute_azimuthal_component(0, nlat, nlon, m, wb, iw, wwbin)

                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do i=1, imm1
                                two1 = wo(i, 2*mp1-1, k)*wts(i)
                                two2 = wo(i, 2*mp1-2, k)*wts(i)
                                tve1 = ve(i, 2*mp1-1, k)*wts(i)
                                tve2 = ve(i, 2*mp1-2, k)*wts(i)
                                do np1=mp1, odd_stride, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*two2 &
                                        +wb(i, np1, iw)*tve1
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*two1 &
                                        -wb(i, np1, iw)*tve2
                                end do
                            end do
                        end do
                        if (odd(nlat))then
                            i = imid
                            do k=1, nt
                                do np1=mp1, odd_stride, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)+wb(i, np1, iw)*ve(i, 2*mp1-1, k)*wts(i)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-wb(i, np1, iw)*ve(i, 2*mp1-2, k)*wts(i)
                                end do
                            end do
                        end if
                    end if

                    if (mp2 > even_stride) cycle main_loop_case_3

                    do k=1, nt
                        do i=1, imm1
                            two1 = wo(i, 2*mp1-1, k)*wts(i)
                            two2 = wo(i, 2*mp1-2, k)*wts(i)
                            tve1 = ve(i, 2*mp1-1, k)*wts(i)
                            tve2 = ve(i, 2*mp1-2, k)*wts(i)
                            do np1=mp2, even_stride, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*tve2 &
                                    +wb(i, np1, iw)*two1
                                bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*tve1 &
                                    -wb(i, np1, iw)*two2
                            end do
                        end do
                    end do

                    if (even(nlat)) cycle main_loop_case_3

                    i = imid

                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*ve(i, 2*mp1-2, k)*wts(i)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*ve(i, 2*mp1-1, k)*wts(i)
                        end do
                    end do
                end do main_loop_case_3
            case(4)
                !
                ! case ityp=4 ,  v even, w odd, and cr and ci equal 0.
                !
                call util%compute_polar_component(1, nlat, nlon, 0, vb, iv, wvbin)
                !
                ! case m=0
                !
                do k=1, nt
                    do i=1, imid
                        tv = ve(i, 1, k)*wts(i)
                        do np1=2, even_stride, 2
                            br(1, np1, k) = br(1, np1, k)+vb(i, np1, iv)*tv
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

                    call util%compute_polar_component(1, nlat, nlon, m, vb, iv, wvbin)
                    call util%compute_azimuthal_component(1, nlat, nlon, m, wb, iw, wwbin)

                    if (mp2 > even_stride) cycle main_loop_case_4

                    do k=1, nt
                        do i=1, imm1
                            two1 = wo(i, 2*mp1-1, k)*wts(i)
                            two2 = wo(i, 2*mp1-2, k)*wts(i)
                            tve1 = ve(i, 2*mp1-1, k)*wts(i)
                            tve2 = ve(i, 2*mp1-2, k)*wts(i)
                            do np1=mp2, even_stride, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*tve2 &
                                    +wb(i, np1, iw)*two1
                                bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*tve1 &
                                    -wb(i, np1, iw)*two2
                            end do
                        end do
                    end do

                    if (even(nlat)) cycle main_loop_case_4

                    i = imid

                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*ve(i, 2*mp1-2, k)*wts(i)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*ve(i, 2*mp1-1, k)*wts(i)
                        end do
                    end do
                end do main_loop_case_4
            case(5)
                !
                ! case ityp=5   v even, w odd, and br and bi equal zero
                !
                call util%compute_polar_component(2, nlat, nlon, 0, vb, iv, wvbin)
                !
                ! case m=0
                !
                do k=1, nt
                    do i=1, imm1
                        tw = wo(i, 1, k)*wts(i)
                        do np1=3, odd_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-vb(i, np1, iv)*tw
                        end do
                    end do
                end do

                ! case m = 1 through nlat-1
                if (mmax < 2) return

                main_loop_case_5: do mp1=2, mmax

                    m = mp1-1
                    mp2 = mp1+1

                    call util%compute_polar_component(2, nlat, nlon, m, vb, iv, wvbin)
                    call util%compute_azimuthal_component(2, nlat, nlon, m, wb, iw, wwbin)

                    if (mp1 > odd_stride) cycle main_loop_case_5

                    do k=1, nt
                        do i=1, imm1
                            two1 = wo(i, 2*mp1-1, k)*wts(i)
                            two2 = wo(i, 2*mp1-2, k)*wts(i)
                            tve1 = ve(i, 2*mp1-1, k)*wts(i)
                            tve2 = ve(i, 2*mp1-2, k)*wts(i)
                            do np1=mp1, odd_stride, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*two2 &
                                    +wb(i, np1, iw)*tve1
                                ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*two1 &
                                    -wb(i, np1, iw)*tve2
                            end do
                        end do
                    end do

                    if (even(nlat)) cycle main_loop_case_5

                    i = imid

                    do k=1, nt
                        do np1=mp1, odd_stride, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)+wb(i, np1, iw)*ve(i, 2*mp1-1, k)*wts(i)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-wb(i, np1, iw)*ve(i, 2*mp1-2, k)*wts(i)
                        end do
                    end do
                end do main_loop_case_5
            case(6)
                !
                ! case ityp=6 ,  v odd , w even
                !
                call util%compute_polar_component(0, nlat, nlon, 0, vb, iv, wvbin)
                !
                ! case m=0
                !
                do k=1, nt
                    do i=1, imid
                        tw = we(i, 1, k)*wts(i)
                        do np1=2, even_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-vb(i, np1, iv)*tw
                        end do
                    end do
                end do

                do k=1, nt
                    do i=1, imm1
                        tv = vo(i, 1, k)*wts(i)
                        do np1=3, odd_stride, 2
                            br(1, np1, k) = br(1, np1, k)+vb(i, np1, iv)*tv
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
                    call util%compute_polar_component(0, nlat, nlon, m, vb, iv, wvbin)
                    call util%compute_azimuthal_component(0, nlat, nlon, m, wb, iw, wwbin)

                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do i=1, imm1
                                twe1 = we(i, 2*mp1-1, k)*wts(i)
                                twe2 = we(i, 2*mp1-2, k)*wts(i)
                                tvo1 = vo(i, 2*mp1-1, k)*wts(i)
                                tvo2 = vo(i, 2*mp1-2, k)*wts(i)
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*tvo2 &
                                        +wb(i, np1, iw)*twe1
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*tvo1 &
                                        -wb(i, np1, iw)*twe2
                                end do
                            end do
                        end do
                        if (odd(nlat)) then
                            i = imid
                            do k=1, nt
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+wb(i, np1, iw)*we(i, 2*mp1-1, k)*wts(i)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)-wb(i, np1, iw)*we(i, 2*mp1-2, k)*wts(i)
                                end do
                            end do
                        end if
                    end if

                    if (mp2 > even_stride) cycle main_loop_case_6

                    do k=1, nt
                        do i=1, imm1
                            twe1 = we(i, 2*mp1-1, k)*wts(i)
                            twe2 = we(i, 2*mp1-2, k)*wts(i)
                            tvo1 = vo(i, 2*mp1-1, k)*wts(i)
                            tvo2 = vo(i, 2*mp1-2, k)*wts(i)
                            do np1=mp2, even_stride, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*twe2 &
                                    +wb(i, np1, iw)*tvo1
                                ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*twe1 &
                                    -wb(i, np1, iw)*tvo2
                            end do
                        end do
                    end do

                    if (even(nlat)) cycle main_loop_case_6

                    i = imid

                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*we(i, 2*mp1-2, k)*wts(i)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*we(i, 2*mp1-1, k)*wts(i)
                        end do
                    end do
                end do main_loop_case_6
            case(7)

                ! case ityp=7   v odd, w even, and cr and ci equal zero
                call util%compute_polar_component(2, nlat, nlon, 0, vb, iv, wvbin)

                ! case m=0
                do k=1, nt
                    do i=1, imm1
                        tv = vo(i, 1, k)*wts(i)
                        do np1=3, odd_stride, 2
                            br(1, np1, k) = br(1, np1, k)+vb(i, np1, iv)*tv
                        end do
                    end do
                end do

                ! case m = 1 through nlat-1
                if (mmax < 2) return

                main_loop_case_7: do mp1=2, mmax

                    m = mp1-1
                    mp2 = mp1+1

                    call util%compute_polar_component(2, nlat, nlon, m, vb, iv, wvbin)
                    call util%compute_azimuthal_component(2, nlat, nlon, m, wb, iw, wwbin)

                    if (mp1 > odd_stride) cycle main_loop_case_7

                    do k=1, nt
                        do i=1, imm1
                            twe1 = we(i, 2*mp1-1, k)*wts(i)
                            twe2 = we(i, 2*mp1-2, k)*wts(i)
                            tvo1 = vo(i, 2*mp1-1, k)*wts(i)
                            tvo2 = vo(i, 2*mp1-2, k)*wts(i)
                            do np1=mp1, odd_stride, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*tvo2 &
                                    +wb(i, np1, iw)*twe1
                                bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*tvo1 &
                                    -wb(i, np1, iw)*twe2
                            end do
                        end do
                    end do

                    if (even(nlat)) cycle main_loop_case_7

                    i = imid

                    do k=1, nt
                        do np1=mp1, odd_stride, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+wb(i, np1, iw)*we(i, 2*mp1-1, k)*wts(i)
                            bi(mp1, np1, k) = bi(mp1, np1, k)-wb(i, np1, iw)*we(i, 2*mp1-2, k)*wts(i)
                        end do
                    end do
                end do main_loop_case_7
            case(8)

                ! case ityp=8   v odd, w even, and both br and bi equal zero
                call util%compute_polar_component(1, nlat, nlon, 0, vb, iv, wvbin)

                ! case m=0
                do k=1, nt
                    do i=1, imid
                        tw = we(i, 1, k)*wts(i)
                        do np1=2, even_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-vb(i, np1, iv)*tw
                        end do
                    end do
                end do

                ! case m = 1 through nlat-1
                if (mmax < 2) return

                main_loop_case_8: do mp1=2, mmax

                    m = mp1-1
                    mp2 = mp1+1

                    call util%compute_polar_component(1, nlat, nlon, m, vb, iv, wvbin)
                    call util%compute_azimuthal_component(1, nlat, nlon, m, wb, iw, wwbin)

                    if (mp2 > even_stride) cycle main_loop_case_8

                    do k=1, nt
                        do i=1, imm1
                            twe1 = we(i, 2*mp1-1, k)*wts(i)
                            twe2 = we(i, 2*mp1-2, k)*wts(i)
                            tvo1 = vo(i, 2*mp1-1, k)*wts(i)
                            tvo2 = vo(i, 2*mp1-2, k)*wts(i)
                            do np1=mp2, even_stride, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*twe2 &
                                    +wb(i, np1, iw)*tvo1
                                ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*twe1 &
                                    -wb(i, np1, iw)*tvo2
                            end do
                        end do
                    end do

                    if (even(nlat)) cycle main_loop_case_8

                    do k=1, nt
                        do  np1=mp2, even_stride, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k) &
                                - vb(imid, np1, iv) * we(imid, 2*mp1-2, k) * wts(imid)
                            ci(mp1, np1, k) = ci(mp1, np1, k) &
                                - vb(imid, np1, iv) * we(imid, 2*mp1-1, k) * wts(imid)
                        end do
                    end do
                end do main_loop_case_8
        end select vector_symmetry_cases

    end subroutine vhagc_lower_utility_routine

end submodule vector_analysis_gaussian_grid
