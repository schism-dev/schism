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
submodule(vector_analysis_routines) vector_analysis_regular_grid_saved

contains

    ! Purpose:
    !
    !     subroutine vhaes(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
    !                      mdab, ndab, wvhaes, lvhaes, work, lwork, ierror)
    !
    !     subroutine vhaes performs the vector spherical harmonic analysis
    !     on the vector field (v, w) and stores the result in the arrays
    !     br, bi, cr, and ci. v(i, j) and w(i, j) are the colatitudinal
    !     (measured from the north pole) and east longitudinal components
    !     respectively, located at colatitude theta(i) = (i-1)*pi/(nlat-1)
    !     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral
    !     representation of (v, w) is given at output parameters v, w in
    !     subroutine vhses.
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
    !     nt     the number of analyses.  in the program that calls vhaes,
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
    !            the program that calls vhaes. if ityp <= 2 then idvw
    !            must be at least nlat.  if ityp > 2 and nlat is
    !            even then idvw must be at least nlat/2. if ityp > 2
    !            and nlat is odd then idvw must be at least (nlat + 1)/2.
    !
    !     jdvw   the second dimension of the arrays v, w as it appears in
    !            the program that calls vhaes. jdvw must be at least nlon.
    !
    !     mdab   the first dimension of the arrays br, bi, cr, and ci as it
    !            appears in the program that calls vhaes. mdab must be at
    !            least min(nlat, nlon/2) if nlon is even or at least
    !            min(nlat, (nlon + 1)/2) if nlon is odd.
    !
    !     ndab   the second dimension of the arrays br, bi, cr, and ci as it
    !            appears in the program that calls vhaes. ndab must be at
    !            least nlat.
    !
    !     lvhaes an array which must be initialized by subroutine vhaesi.
    !            once initialized, wvhaes can be used repeatedly by vhaes
    !            as long as nlon and nlat remain unchanged.  wvhaes must
    !            not be altered between calls of vhaes.
    !
    !     lvhaes the dimension of the array wvhaes as it appears in the
    !            program that calls vhaes. define
    !
    !               l1 = min(nlat, nlon/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lvhaes must be at least
    !
    !            l1*l2(2*nlat-l1+1)+nlon+15
    !
    !     output parameters
    !
    !     br, bi  two or three dimensional arrays (see input parameter nt)
    !     cr, ci  that contain the vector spherical harmonic coefficients
    !            in the spectral representation of v(i, j) and w(i, j) given
    !            in the discription of subroutine vhses. br(mp1, np1),
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
    !            = 9  error in the specification of lvhaes
    !
    module subroutine vhaes(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
        mdab, ndab, wvhaes, ierror)

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
        real(wp),    intent(in)  :: wvhaes(:)
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip) :: idv, imid, idz, ist, lnl, lzimn, mmax
        integer(ip) :: required_wavetable_size
        integer(ip) :: lwork, workspace_indices(6)
        type(VectorAnalysisUtility) :: util

        imid = (nlat + 1)/2
        mmax = min(nlat, (nlon + 1)/2)
        idz = (mmax*(2*nlat-mmax+1))/2
        lzimn = idz*imid
        required_wavetable_size = (2 * lzimn) + nlon + 15

        ! Check calling arguments
        call util%check_vector_analysis_inputs(nlat, nlon, ityp, idvw, jdvw, &
            mdab, ndab, nt, required_wavetable_size, wvhaes, ierror)

        ! Check error flag
        if (ierror /= 0) return

        select case (ityp)
            case(0:2)
                idv = nlat
                ist = imid
            case default
                idv = imid
                ist = 0
        end select

        lnl = nt*idv*nlon

        ! Set required workspace size
        lwork = 2*lnl+idv*nlon

        block
            real(wp) :: work(lwork)

            !  Set workspace indices
            workspace_indices = get_vhaes_workspace_indices(ist, lnl, lzimn)

            associate (&
                iw1 => workspace_indices(1), &
                iw2 => workspace_indices(2), &
                iw3 => workspace_indices(3), &
                iw4 => workspace_indices(4), &
                jw1 => workspace_indices(5), &
                jw2 => workspace_indices(6) &
                )
                call vhaes_lower_utility_routine(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, ndab, &
                    br, bi, cr, ci, idv, work, work(iw1:), work(iw2:), work(iw3:), &
                    work(iw4:), idz, wvhaes, wvhaes(jw1:), wvhaes(jw2:))
            end associate
        end block

    end subroutine vhaes

    ! Purpose:
    !
    !     subroutine vhaesi(nlat, nlon, wvhaes, ierror)
    !
    !     subroutine vhaesi initializes the array wvhaes which can then be
    !     used repeatedly by subroutine vhaes until nlat or nlon is changed.
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
    !     lvhaes the dimension of the array wvhaes as it appears in the
    !            program that calls vhaes. define
    !
    !               l1 = min(nlat, nlon/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lvhaes must be at least
    !
    !               l1*l2*(2*nlat-l1+1)+nlon+15
    !
    !     output parameters
    !
    !     wvhaes an array which is initialized for use by subroutine vhaes.
    !            once initialized, wvhaes can be used repeatedly by vhaes
    !            as long as nlat or nlon remain unchanged.  wvhaes must not
    !            be altered between calls of vhaes.
    !
    !
    !     ierror = 0  no errors
    !            = 1  error in the specification of nlat
    !            = 2  error in the specification of nlon
    !            = 3  error in the specification of lvhaes
    !
    module subroutine vhaesi(nlat, nlon, wvhaes, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(out) :: wvhaes(:)
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip) :: imid, labc, lzimn, mmax
        integer(ip) :: required_wavetable_size
        integer(ip) :: workspace_indices(4)
        integer(ip) :: lwork, ldwork
        type(SpherepackUtility) :: util

        mmax = min(nlat, (nlon + 1)/2)
        imid = (nlat + 1)/2
        lzimn = (imid*mmax*(2*nlat-mmax+1))/2
        labc = 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
        required_wavetable_size = 2*lzimn+nlon+15

        ! Check calling arguments
        if (nlat < 3) then
            ierror = 1
        else if (nlon < 1) then
            ierror = 2
        else if (size(wvhaes) < required_wavetable_size) then
            ierror = 3
        else
            ierror = 0
        end if

        ! Check error flag
        if (ierror /= 0) return

        !  Set workspace indices
        workspace_indices = get_vhaesi_workspace_indices(lzimn, nlat, imid, mmax)

        ! Set required workspace sizes
        lwork = 5*nlat*imid+labc
        ldwork = 2 * (nlat + 1)

        block
            real(wp) :: work(lwork), dwork(ldwork)

            associate (&
                jw1 => workspace_indices(1), &
                jw2 => workspace_indices(2), &
                iw1 => workspace_indices(3), &
                idz => workspace_indices(4) &
                )
                call vhaesi_lower_utility_routine( &
                    nlat, nlon, imid, wvhaes, wvhaes(jw1:), idz, work, work(iw1:), dwork)

                call util%hfft%initialize(nlon, wvhaes(jw2:))
            end associate
        end block

    end subroutine vhaesi

    pure function get_vhaes_workspace_indices(ist, lnl, lzimn) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: ist
        integer(ip), intent(in)  :: lnl
        integer(ip), intent(in)  :: lzimn
        integer(ip)              :: return_value(6)

        associate (i => return_value)
            i(1) = ist+1
            i(2) = lnl+1
            i(3) = i(2)+ist
            i(4) = i(2)+lnl
            i(5) = lzimn+1
            i(6) = i(5)+lzimn
        end associate

    end function get_vhaes_workspace_indices

    subroutine vhaes_lower_utility_routine(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
        ndab, br, bi, cr, ci, idv, ve, vo, we, wo, work, idz, zv, zw, wrfft)

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
        real(wp),    intent(out) :: work(*)
        integer(ip), intent(in)  :: idz
        real(wp),    intent(in)  :: zv(idz, *)
        real(wp),    intent(in)  :: zw(idz, *)
        real(wp),    intent(in)  :: wrfft(:)

        ! Local variables
        integer(ip) :: i, imm1, k, m, mb
        integer(ip) :: mmax, mp1, mp2, odd_stride, even_stride, np1
        type(VectorAnalysisUtility) :: util

        call util%analysis_setup(idvw, jdvw, mdab, ndab, &
            bi, br, ci, cr, even_stride, idv, imid, imm1, ityp, &
            mmax, nlat, nlon, nt, odd_stride,  v, ve, vo, w, we, wo, wrfft)

        !  Compute coefficients br, bi, cr, ci
        select case (ityp)
            case (0)

                !  case ityp=0,  no symmetries

                ! case m=0
                do k=1, nt
                    do i=1, imid
                        do np1=2, even_stride, 2
                            br(1, np1, k) = br(1, np1, k)+zv(np1, i)*ve(i, 1, k)
                            cr(1, np1, k) = cr(1, np1, k)-zv(np1, i)*we(i, 1, k)
                        end do
                    end do
                end do

                do k=1, nt
                    do i=1, imm1
                        do np1=3, odd_stride, 2
                            br(1, np1, k) = br(1, np1, k)+zv(np1, i)*vo(i, 1, k)
                            cr(1, np1, k) = cr(1, np1, k)-zv(np1, i)*wo(i, 1, k)
                        end do
                    end do
                end do

                ! case m = 1 through nlat-1
                if (mmax < 2) return

                main_loop_case_0: do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1

                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do i=1, imm1
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, i)*vo(i, 2*mp1-2, k) &
                                        +zw(np1+mb, i)*we(i, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, i)*vo(i, 2*mp1-1, k) &
                                        -zw(np1+mb, i)*we(i, 2*mp1-2, k)
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, i)*wo(i, 2*mp1-2, k) &
                                        +zw(np1+mb, i)*ve(i, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, i)*wo(i, 2*mp1-1, k) &
                                        -zw(np1+mb, i)*ve(i, 2*mp1-2, k)
                                end do
                            end do
                        end do
                        if (mod(nlat, 2) /= 0) then
                            do k=1, nt
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+zw(np1+mb, imid)*we(imid, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)-zw(np1+mb, imid)*we(imid, 2*mp1-2, k)
                                    cr(mp1, np1, k) = cr(mp1, np1, k)+zw(np1+mb, imid)*ve(imid, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-zw(np1+mb, imid)*ve(imid, 2*mp1-2, k)
                                end do
                            end do
                        end if
                    end if

                    if (mp2 > even_stride) cycle main_loop_case_0

                    do k=1, nt
                        do i=1, imm1
                            do np1=mp2, even_stride, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, i)*ve(i, 2*mp1-2, k) &
                                    +zw(np1+mb, i)*wo(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, i)*ve(i, 2*mp1-1, k) &
                                    -zw(np1+mb, i)*wo(i, 2*mp1-2, k)
                                cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, i)*we(i, 2*mp1-2, k) &
                                    +zw(np1+mb, i)*vo(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, i)*we(i, 2*mp1-1, k) &
                                    -zw(np1+mb, i)*vo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do

                    if (mod(nlat, 2) == 0) cycle main_loop_case_0

                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, imid)*ve(imid, 2*mp1-2, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, imid)*ve(imid, 2*mp1-1, k)
                            cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, imid)*we(imid, 2*mp1-2, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, imid)*we(imid, 2*mp1-1, k)
                        end do
                    end do
                end do main_loop_case_0
            case (1)
                !
                !  case ityp=1 ,  no symmetries but cr and ci equal zero
                !
                !    case m=0
                !
                do k=1, nt
                    do i=1, imid
                        do np1=2, even_stride, 2
                            br(1, np1, k) = br(1, np1, k)+zv(np1, i)*ve(i, 1, k)
                        end do
                    end do
                end do

                do k=1, nt
                    do i=1, imm1
                        do np1=3, odd_stride, 2
                            br(1, np1, k) = br(1, np1, k)+zv(np1, i)*vo(i, 1, k)
                        end do
                    end do
                end do
                !
                !  case m = 1 through nlat-1
                !
                if (mmax < 2) return

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1

                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do i=1, imm1
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, i)*vo(i, 2*mp1-2, k) &
                                        +zw(np1+mb, i)*we(i, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, i)*vo(i, 2*mp1-1, k) &
                                        -zw(np1+mb, i)*we(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mod(nlat, 2) /= 0) then
                            do k=1, nt
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+zw(np1+mb, imid)*we(imid, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)-zw(np1+mb, imid)*we(imid, 2*mp1-2, k)
                                end do
                            end do
                        end if
                    end if

                    if (mp2 > even_stride) return

                    do k=1, nt
                        do i=1, imm1
                            do np1=mp2, even_stride, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, i)*ve(i, 2*mp1-2, k) &
                                    +zw(np1+mb, i)*wo(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, i)*ve(i, 2*mp1-1, k) &
                                    -zw(np1+mb, i)*wo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do

                    if (mod(nlat, 2) == 0) return

                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, imid)*ve(imid, 2*mp1-2, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, imid)*ve(imid, 2*mp1-1, k)
                        end do
                    end do
                end do
            case (2)
                !
                !  case ityp=2 ,  no symmetries but br and bi equal zero
                !
                !    case m=0
                !
                do k=1, nt
                    do i=1, imid
                        do np1=2, even_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-zv(np1, i)*we(i, 1, k)
                        end do
                    end do
                end do

                do k=1, nt
                    do i=1, imm1
                        do np1=3, odd_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-zv(np1, i)*wo(i, 1, k)
                        end do
                    end do
                end do
                !
                !  case m = 1 through nlat-1
                !
                if (mmax < 2) return

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1

                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do i=1, imm1
                                do np1=mp1, odd_stride, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, i)*wo(i, 2*mp1-2, k) &
                                        +zw(np1+mb, i)*ve(i, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, i)*wo(i, 2*mp1-1, k) &
                                        -zw(np1+mb, i)*ve(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mod(nlat, 2) /= 0) then
                            do k=1, nt
                                do np1=mp1, odd_stride, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)+zw(np1+mb, imid)*ve(imid, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-zw(np1+mb, imid)*ve(imid, 2*mp1-2, k)
                                end do
                            end do
                        end if
                    end if

                    if (mp2 > even_stride) return

                    do k=1, nt
                        do i=1, imm1
                            do np1=mp2, even_stride, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, i)*we(i, 2*mp1-2, k) &
                                    +zw(np1+mb, i)*vo(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, i)*we(i, 2*mp1-1, k) &
                                    -zw(np1+mb, i)*vo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do

                    if (mod(nlat, 2) == 0) return

                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, imid)*we(imid, 2*mp1-2, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, imid)*we(imid, 2*mp1-1, k)
                        end do
                    end do
                end do
            case (3)
                !
                !  case ityp=3 ,  v even , w odd
                !
                !    case m=0
                !
                do k=1, nt
                    do i=1, imid
                        do np1=2, even_stride, 2
                            br(1, np1, k) = br(1, np1, k)+zv(np1, i)*ve(i, 1, k)
                        end do
                    end do
                end do

                do k=1, nt
                    do i=1, imm1
                        do np1=3, odd_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-zv(np1, i)*wo(i, 1, k)
                        end do
                    end do
                end do
                !
                !  case m = 1 through nlat-1
                !
                if (mmax < 2) return

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1

                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do i=1, imm1
                                do np1=mp1, odd_stride, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, i)*wo(i, 2*mp1-2, k) &
                                        +zw(np1+mb, i)*ve(i, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, i)*wo(i, 2*mp1-1, k) &
                                        -zw(np1+mb, i)*ve(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mod(nlat, 2) /= 0) then
                            do k=1, nt
                                do np1=mp1, odd_stride, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)+zw(np1+mb, imid)*ve(imid, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-zw(np1+mb, imid)*ve(imid, 2*mp1-2, k)
                                end do
                            end do
                        end if
                    end if

                    if (mp2 > even_stride) return

                    do k=1, nt
                        do i=1, imm1
                            do np1=mp2, even_stride, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, i)*ve(i, 2*mp1-2, k) &
                                    +zw(np1+mb, i)*wo(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, i)*ve(i, 2*mp1-1, k) &
                                    -zw(np1+mb, i)*wo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do

                    if (mod(nlat, 2) == 0) return

                    do  k=1, nt
                        do  np1=mp2, even_stride, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, imid)*ve(imid, 2*mp1-2, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, imid)*ve(imid, 2*mp1-1, k)
                        end do
                    end do
                end do
            case (4)
                !
                !  case ityp=4 ,  v even, w odd, and cr and ci equal 0.
                !
                !    case m=0
                !
                do k=1, nt
                    do i=1, imid
                        do np1=2, even_stride, 2
                            br(1, np1, k) = br(1, np1, k)+zv(np1, i)*ve(i, 1, k)
                        end do
                    end do
                end do
                !
                !  case m = 1 through nlat-1
                !
                if (mmax < 2) return

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1

                    if (mp2 > even_stride) exit

                    do k=1, nt
                        do i=1, imm1
                            do np1=mp2, even_stride, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, i)*ve(i, 2*mp1-2, k) &
                                    +zw(np1+mb, i)*wo(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, i)*ve(i, 2*mp1-1, k) &
                                    -zw(np1+mb, i)*wo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do

                    if (mod(nlat, 2) == 0) return

                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, imid)*ve(imid, 2*mp1-2, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, imid)*ve(imid, 2*mp1-1, k)
                        end do
                    end do
                end do
            case (5)
                !
                !  case ityp=5   v even, w odd, and br and bi equal zero
                !
                !    case m=0
                !
                do k=1, nt
                    do i=1, imm1
                        do np1=3, odd_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-zv(np1, i)*wo(i, 1, k)
                        end do
                    end do
                end do
                !
                !  case m = 1 through nlat-1
                !
                if (mmax < 2) return

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1

                    if (mp1 > odd_stride) exit

                    do k=1, nt
                        do i=1, imm1
                            do np1=mp1, odd_stride, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, i)*wo(i, 2*mp1-2, k) &
                                    +zw(np1+mb, i)*ve(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, i)*wo(i, 2*mp1-1, k) &
                                    -zw(np1+mb, i)*ve(i, 2*mp1-2, k)
                            end do
                        end do
                    end do

                    if (mod(nlat, 2) == 0) exit

                    do k=1, nt
                        do np1=mp1, odd_stride, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)+zw(np1+mb, imid)*ve(imid, 2*mp1-1, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-zw(np1+mb, imid)*ve(imid, 2*mp1-2, k)
                        end do
                    end do
                end do
            case (6)
                !
                !  case ityp=6 ,  v odd , w even
                !
                !    case m=0
                !
                do k=1, nt
                    do i=1, imid
                        do np1=2, even_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-zv(np1, i)*we(i, 1, k)
                        end do
                    end do
                end do

                do k=1, nt
                    do i=1, imm1
                        do np1=3, odd_stride, 2
                            br(1, np1, k) = br(1, np1, k)+zv(np1, i)*vo(i, 1, k)
                        end do
                    end do
                end do
                !
                !  case m = 1 through nlat-1
                !
                if (mmax < 2) return

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1

                    if (mp1 <= odd_stride) then

                        do  k=1, nt
                            do  i=1, imm1
                                do  np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, i)*vo(i, 2*mp1-2, k) &
                                        +zw(np1+mb, i)*we(i, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, i)*vo(i, 2*mp1-1, k) &
                                        -zw(np1+mb, i)*we(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mod(nlat, 2) /= 0) then
                            do  k=1, nt
                                do  np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+zw(np1+mb, imid)*we(imid, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)-zw(np1+mb, imid)*we(imid, 2*mp1-2, k)
                                end do
                            end do
                        end if
                    end if

                    if (mp2 > even_stride) return

                    do k=1, nt
                        do  i=1, imm1
                            do  np1=mp2, even_stride, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, i)*we(i, 2*mp1-2, k) &
                                    +zw(np1+mb, i)*vo(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, i)*we(i, 2*mp1-1, k) &
                                    -zw(np1+mb, i)*vo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do

                    if (mod(nlat, 2) == 0) return

                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, imid)*we(imid, 2*mp1-2, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, imid)*we(imid, 2*mp1-1, k)
                        end do
                    end do
                end do
            case (7)
                !
                !  case ityp=7   v odd, w even, and cr and ci equal zero
                !
                !    case m=0
                !
                do k=1, nt
                    do i=1, imm1
                        do np1=3, odd_stride, 2
                            br(1, np1, k) = br(1, np1, k)+zv(np1, i)*vo(i, 1, k)
                        end do
                    end do
                end do
                !
                !  case m = 1 through nlat-1
                !
                if (mmax < 2) return

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1

                    if (mp1 > odd_stride) return

                    do  k=1, nt
                        do  i=1, imm1
                            do  np1=mp1, odd_stride, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, i)*vo(i, 2*mp1-2, k) &
                                    +zw(np1+mb, i)*we(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, i)*vo(i, 2*mp1-1, k) &
                                    -zw(np1+mb, i)*we(i, 2*mp1-2, k)
                            end do
                        end do
                    end do

                    if (mod(nlat, 2) == 0) return

                    do k=1, nt
                        do np1=mp1, odd_stride, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+zw(np1+mb, imid)*we(imid, 2*mp1-1, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)-zw(np1+mb, imid)*we(imid, 2*mp1-2, k)
                        end do
                    end do
                end do
            case (8)
                !
                !  case ityp=8   v odd, w even, and both br and bi equal zero
                !
                !    case m=0
                !
                do k=1, nt
                    do i=1, imid
                        do np1=2, even_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-zv(np1, i)*we(i, 1, k)
                        end do
                    end do
                end do
                !
                ! case m = 1 through nlat-1
                !
                if (mmax < 2) return

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1

                    if (mp2 > even_stride) return

                    do k=1, nt
                        do i=1, imm1
                            do np1=mp2, even_stride, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, i)*we(i, 2*mp1-2, k) &
                                    +zw(np1+mb, i)*vo(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, i)*we(i, 2*mp1-1, k) &
                                    -zw(np1+mb, i)*vo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do

                    if (mod(nlat, 2) == 0) return

                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, imid)*we(imid, 2*mp1-2, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, imid)*we(imid, 2*mp1-1, k)
                        end do
                    end do
                end do
        end select

    end subroutine vhaes_lower_utility_routine

    pure function get_vhaesi_workspace_indices(lzimn, nlat, imid, mmax) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: lzimn
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: imid
        integer(ip), intent(in)  :: mmax
        integer(ip)              :: return_value(4)


        associate (i => return_value)
            i(1) = lzimn+1
            i(2) = 2*lzimn+1
            i(3) = 3*nlat*imid+1
            i(4) = (mmax*(2*nlat-mmax+1))/2
        end associate

    end function get_vhaesi_workspace_indices

    subroutine vhaesi_lower_utility_routine(nlat, nlon, imid, zv, zw, idz, zin, wzvin, dwork)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: imid
        real(wp),    intent(out) :: zv(idz, *)
        real(wp),    intent(out) :: zw(idz, *)
        integer(ip), intent(in)  :: idz
        real(wp),    intent(out) :: zin(imid, nlat, 3)
        real(wp),    intent(out) :: wzvin(*)
        real(wp),    intent(out) :: dwork(*)

        ! Local variables
        integer(ip)         :: i3, m, mn, mp1, np1, mmax
        type(SpherepackUtility) :: util

        mmax = min(nlat, (nlon + 1)/2)

        call util%zvinit(nlat, nlon, wzvin, dwork)

        do mp1=1, mmax
            m = mp1-1
            call util%zvin(0, nlat, nlon, m, zin, i3, wzvin)
            do np1=mp1, nlat
                mn = m*(nlat-1)-(m*(m-1))/2+np1
                zv(mn, 1:imid) = zin(1:imid, np1, i3)
            end do
        end do

        call util%zwinit(nlat, nlon, wzvin, dwork)

        do mp1=1, mmax
            m = mp1-1
            call util%zwin(0, nlat, nlon, m, zin, i3, wzvin)
            do np1=mp1, nlat
                mn = m*(nlat-1)-(m*(m-1))/2+np1
                zw(mn, 1:imid) = zin(1:imid, np1, i3)
            end do
        end do

    end subroutine vhaesi_lower_utility_routine

end submodule vector_analysis_regular_grid_saved
