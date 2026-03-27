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
submodule(vector_analysis_routines) vector_analysis_gaussian_grid_saved

contains

    ! Purpose:
    !
    !     subroutine vhags(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
    !                      mdab, ndab, wvhags, ierror)
    !
    !     subroutine vhags performs the vector spherical harmonic analysis
    !     on the vector field (v, w) and stores the result in the arrays
    !     br, bi, cr, and ci. v(i, j) and w(i, j) are the colatitudinal
    !     (measured from the north pole) and east longitudinal components
    !     respectively, located at the gaussian colatitude point theta(i)
    !     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral
    !     representation of (v, w) is given at output parameters v, w in
    !     subroutine vhses.
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
    !     nt     the number of analyses.  in the program that calls vhags,
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
    !            components at the gaussian colatitude point theta(i)
    !            and longitude phi(j) = (j-1)*2*pi/nlon. the index ranges
    !            are defined above at the input parameter ityp.
    !
    !     idvw   the first dimension of the arrays v, w as it appears in
    !            the program that calls vhags. if ityp <= 2 then idvw
    !            must be at least nlat.  if ityp > 2 and nlat is
    !            even then idvw must be at least nlat/2. if ityp > 2
    !            and nlat is odd then idvw must be at least (nlat + 1)/2.
    !
    !     jdvw   the second dimension of the arrays v, w as it appears in
    !            the program that calls vhags. jdvw must be at least nlon.
    !
    !     mdab   the first dimension of the arrays br, bi, cr, and ci as it
    !            appears in the program that calls vhags. mdab must be at
    !            least min(nlat, nlon/2) if nlon is even or at least
    !            min(nlat, (nlon + 1)/2) if nlon is odd.
    !
    !     ndab   the second dimension of the arrays br, bi, cr, and ci as it
    !            appears in the program that calls vhags. ndab must be at
    !            least nlat.
    !
    !     wvhags an array which must be initialized by subroutine vhgsi.
    !            once initialized, wvhags can be used repeatedly by vhags
    !            as long as nlon and nlat remain unchanged.  wvhags must
    !            not be altered between calls of vhags.
    !
    !     lvhags the dimension of the array wvhags as it appears in the
    !            program that calls vhags. define
    !
    !               l1 = min(nlat, nlon/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lvhags must be at least
    !
    !            l1*l2(nlat+ nlat-l1+1)+ nlon + 15
    !
    !        ??? (nlat + 1)*(nlat + 1)*nlat/2 + nlon + 15
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
    !            = 9  error in the specification of lvhags
    !
    module subroutine vhags(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
        mdab, ndab, wvhags, ierror)

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
        real(wp),    intent(in)  :: wvhags(:)
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip) :: idv, idz, imid, ist
        integer(ip) :: required_wavetable_size
        integer(ip) :: lnl, lzimn, mmax, lwork
        integer(ip) :: workspace_indices(7)
        type(VectorAnalysisUtility) :: util

        mmax = min(nlat, (nlon + 1)/2)
        idz = (mmax * ((2 * nlat) - mmax + 1))/2
        imid = (nlat + 1)/2
        lzimn = idz * imid
        required_wavetable_size = (2 * lzimn) + nlon + 15

        ! Check calling arguments
        call util%check_vector_analysis_inputs(nlat, nlon, ityp, idvw, jdvw, &
            mdab, ndab, nt, required_wavetable_size, wvhags, ierror)

        ! Check error flag
        if (ierror /= 0) return

        select case(ityp)
            case(0:2)
                ist = imid
                idv = nlat
            case default
                ist = 0
                idv = imid
        end select

        lnl = nt*idv*nlon

        ! Set required workspace size
        lwork = 2*lnl+idv*nlon

        block
            real(wp) :: work(lwork)

            ! Compute workspace pointers
            workspace_indices = get_vhags_workspace_indices(nlat, imid, ist, lnl)

            associate (&
                jw1 => workspace_indices(1), &
                jw2 => workspace_indices(2), &
                jw3 => workspace_indices(3), &
                iw1 => workspace_indices(4), &
                iw2 => workspace_indices(5), &
                iw3 => workspace_indices(6), &
                iw4 => workspace_indices(7) &
                )
                call vhags_lower_utility_routine(nlat, nlon, ityp, nt, imid, idvw, jdvw, &
                    v, w, mdab, ndab, br, bi, cr, ci, idv, work, work(iw1:), work(iw2:), &
                    work(iw3:), work(iw4:), idz, wvhags(jw1:), wvhags(jw2:), wvhags(jw3:))
            end associate
        end block

    end subroutine vhags

    ! Purpose:
    !
    !     subroutine vhagsi(nlat, nlon, wvhags, ierror)
    !
    !     subroutine vhagsi initializes the array wvhags which can then be
    !     used repeatedly by subroutine vhags until nlat or nlon is changed.
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
    !     lvhags the dimension of the array wvhags as it appears in the
    !            program that calls vhagsi. lvhags must be at least
    !
    !               3*nlat*(nlat + 1)+2  (required by vhagsi)
    !
    !     output parameters
    !
    !     wvhags an array which is initialized for use by subroutine vhags.
    !            once initialized, wvhags can be used repeatedly by vhags
    !            as long as nlat and nlon remain unchanged.  wvhags must not
    !            be altered between calls of vhags.
    !
    !
    !     ierror = 0  no errors
    !            = 1  error in the specification of nlat
    !            = 2  error in the specification of nlon
    !            = 3  error in the specification of lvhags
    !
    module subroutine vhagsi(nlat, nlon, wvhags, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(out) :: wvhags(:)
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip)  :: imid, lmn, ldwork, required_wavetable_size
        integer(ip)  :: workspace_indices(7)
        type(SpherepackUtility) :: util


        imid = (nlat + 1)/2
        lmn = (nlat*(nlat + 1))/2
        required_wavetable_size = 2*(imid*lmn)+ nlon + 15

        !  Check calling arguments
        if (nlat < 3) then
            ierror = 1
        else if (nlon < 1) then
            ierror = 2
        else if (size(wvhags) < required_wavetable_size) then
            ierror = 3
        else
            ierror = 0
        end if

        ! Check error flag
        if (ierror /= 0) return

        !  Compute workspace indices
        workspace_indices = get_vhagsi_workspace_indices(nlat, imid, lmn)

        associate (&
            jw1 => workspace_indices(1), &
            jw2 => workspace_indices(2), &
            jw3 => workspace_indices(3), &
            iw1 => workspace_indices(4), &
            iw2 => workspace_indices(5), &
            iw3 => workspace_indices(6), &
            iw4 => workspace_indices(7) &
            )

            ! Set required workspace size
            ldwork = (nlat*(3*nlat+9)+2)/2
            block
                real(wp) :: dwork(ldwork)

                call precompute_associated_legendre_functions( &
                    nlat, imid, wvhags(jw1:), wvhags(jw2:), &
                    dwork(iw1:), dwork(iw2:), dwork(iw3:), dwork(iw4:))

                call util%hfft%initialize(nlon, wvhags(jw3:))
            end block
        end associate

    end subroutine vhagsi

    pure function get_vhags_workspace_indices(nlat, imid, ist, lnl) &
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
            !
            !  set wvhags pointers
            !
            lmn = nlat*(nlat + 1)/2
            i(1) = 1
            i(2) = i(1)+imid*lmn
            i(3) = i(2)+imid*lmn
            !
            !  set work pointers
            !
            i(4) = ist+1
            i(5) = lnl+1
            i(6) = i(5)+ist
            i(7) = i(5)+lnl
        end associate

    end function get_vhags_workspace_indices

    subroutine vhags_lower_utility_routine(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
        ndab, br, bi, cr, ci, idv, ve, vo, we, wo, work, idz, vb, wb, wrfft)

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
        real(wp),    intent(in)  :: vb(imid, *)
        real(wp),    intent(in)  :: wb(imid, *)
        real(wp),    intent(in)  :: wrfft(:)

        ! Local variables
        integer(ip) :: i, imm1, k, m, mb, mmax, mp1, mp2
        integer(ip) :: odd_stride, even_stride, np1
        type(VectorAnalysisUtility) :: util

        call util%analysis_setup(idvw, jdvw, mdab, ndab, &
            bi, br, ci, cr, even_stride, idv, imid, imm1, ityp, &
            mmax, nlat, nlon, nt, odd_stride,  v, ve, vo, w, we, wo, wrfft)

        vector_symmetry_cases: select case (ityp)
            case (0)

                ! case ityp = 0,  no symmetries

                ! case m=0
                do k=1, nt
                    do i=1, imid
                        do np1=2, even_stride, 2
                            br(1, np1, k) = br(1, np1, k)+vb(i, np1)*ve(i, 1, k)
                            cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*we(i, 1, k)
                        end do
                    end do
                end do

                do k=1, nt
                    do i=1, imm1
                        do np1=3, odd_stride, 2
                            br(1, np1, k) = br(1, np1, k)+vb(i, np1)*vo(i, 1, k)
                            cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*wo(i, 1, k)
                        end do
                    end do
                end do
                !
                !  case m = 1 through nlat-1
                !
                if (mmax < 2) return

                main_loop_case_0: do mp1=2, mmax
                    m = mp1-1
                    mb = m*nlat-(m*(m + 1))/2
                    mp2 = mp1+1

                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do i=1, imm1
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*we(i, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*we(i, 2*mp1-2, k)
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*ve(i, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*ve(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mod(nlat, 2) /= 0) then
                            do k=1, nt
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+wb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)-wb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                                    cr(mp1, np1, k) = cr(mp1, np1, k)+wb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-wb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                                end do
                            end do
                        end if
                    end if

                    if (mp2 > even_stride) cycle main_loop_case_0

                    do k=1, nt
                        do i=1, imm1
                            do np1=mp2, even_stride, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-2, k) &
                                    +wb(i, np1+mb)*wo(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-1, k) &
                                    -wb(i, np1+mb)*wo(i, 2*mp1-2, k)
                                cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-2, k) &
                                    +wb(i, np1+mb)*vo(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-1, k) &
                                    -wb(i, np1+mb)*vo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do

                    if (mod(nlat, 2) == 0) cycle main_loop_case_0

                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                            cr(mp1, np1, k) = cr(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-1, k)
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
                            br(1, np1, k) = br(1, np1, k)+vb(i, np1)*ve(i, 1, k)
                        end do
                    end do
                end do

                do k=1, nt
                    do i=1, imm1
                        do np1=3, odd_stride, 2
                            br(1, np1, k) = br(1, np1, k)+vb(i, np1)*vo(i, 1, k)
                        end do
                    end do
                end do
                !
                !  case m = 1 through nlat-1
                !
                if (mmax < 2) return

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*nlat-(m*(m + 1))/2
                    mp2 = mp1+1

                    if (mp1 <= odd_stride) then

                        do k=1, nt
                            do i=1, imm1
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*we(i, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*we(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mod(nlat, 2) /= 0) then
                            do k=1, nt
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+wb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)-wb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                                end do
                            end do
                        end if
                    end if

                    if (mp2 > even_stride) return

                    do k=1, nt
                        do i=1, imm1
                            do np1=mp2, even_stride, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-2, k) &
                                    +wb(i, np1+mb)*wo(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-1, k) &
                                    -wb(i, np1+mb)*wo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do

                    if (mod(nlat, 2) == 0) return

                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
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
                            cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*we(i, 1, k)
                        end do
                    end do
                end do

                do k=1, nt
                    do i=1, imm1
                        do np1=3, odd_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*wo(i, 1, k)
                        end do
                    end do
                end do
                !
                !  case m = 1 through nlat-1
                !
                if (mmax < 2) return

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*nlat-(m*(m + 1))/2
                    mp2 = mp1+1

                    if (mp1 <= odd_stride) then

                        do k=1, nt
                            do i=1, imm1
                                do np1=mp1, odd_stride, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*ve(i, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*ve(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mod(nlat, 2) /= 0) then
                            do k=1, nt
                                do np1=mp1, odd_stride, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)+wb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-wb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                                end do
                            end do
                        end if
                    end if

                    if (mp2 > even_stride) return

                    do k=1, nt
                        do i=1, imm1
                            do np1=mp2, even_stride, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-2, k) &
                                    +wb(i, np1+mb)*vo(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-1, k) &
                                    -wb(i, np1+mb)*vo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do

                    if (mod(nlat, 2) == 0) return

                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-1, k)
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
                            br(1, np1, k) = br(1, np1, k)+vb(i, np1)*ve(i, 1, k)
                        end do
                    end do
                end do

                do k=1, nt
                    do i=1, imm1
                        do np1=3, odd_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*wo(i, 1, k)
                        end do
                    end do
                end do
                !
                !  case m = 1 through nlat-1
                !
                if (mmax < 2) return

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*nlat-(m*(m + 1))/2
                    mp2 = mp1+1

                    if (mp1 <= odd_stride) then

                        do k=1, nt
                            do i=1, imm1
                                do np1=mp1, odd_stride, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*ve(i, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*ve(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mod(nlat, 2) /= 0) then

                            do k=1, nt
                                do np1=mp1, odd_stride, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)+wb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-wb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                                end do
                            end do
                        end if
                    end if

                    if (mp2 > even_stride) return

                    do k=1, nt
                        do i=1, imm1
                            do np1=mp2, even_stride, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-2, k) &
                                    +wb(i, np1+mb)*wo(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-1, k) &
                                    -wb(i, np1+mb)*wo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do

                    if (mod(nlat, 2) == 0) return

                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                        end do
                    end do
                end do
            case (4)
                !
                !  case ityp=4 ,  v even, w odd, and cr and ci equal 0.
                !
                ! case m=0
                !
                do k=1, nt
                    do i=1, imid
                        do np1=2, even_stride, 2
                            br(1, np1, k) = br(1, np1, k)+vb(i, np1)*ve(i, 1, k)
                        end do
                    end do
                end do
                !
                !  case m = 1 through nlat-1
                !
                if (mmax < 2) return

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*nlat-(m*(m + 1))/2
                    mp2 = mp1+1

                    if (mp2 > even_stride) return

                    do k=1, nt
                        do i=1, imm1
                            do np1=mp2, even_stride, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-2, k) &
                                    +wb(i, np1+mb)*wo(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-1, k) &
                                    -wb(i, np1+mb)*wo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do

                    if (mod(nlat, 2) == 0) return

                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                        end do
                    end do
                end do
            case (5)
                !
                !  case ityp=5   v even, w odd, and br and bi equal zero
                !
                ! case m=0
                !
                do k=1, nt
                    do i=1, imm1
                        do np1=3, odd_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*wo(i, 1, k)
                        end do
                    end do
                end do
                !
                !  case m = 1 through nlat-1
                !
                if (mmax < 2) return

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*nlat-(m*(m + 1))/2
                    mp2 = mp1+1

                    if (mp1 > odd_stride) return

                    do  k=1, nt
                        do  i=1, imm1
                            do  np1=mp1, odd_stride, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-2, k) &
                                    +wb(i, np1+mb)*ve(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-1, k) &
                                    -wb(i, np1+mb)*ve(i, 2*mp1-2, k)
                            end do
                        end do
                    end do

                    if (mod(nlat, 2) == 0) return

                    do k=1, nt
                        do np1=mp1, odd_stride, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)+wb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-wb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                        end do
                    end do
                end do
            case (6)
                !
                !  case ityp=6 ,  v odd , w even
                !
                ! case m=0
                !
                do k=1, nt
                    do i=1, imid
                        do np1=2, even_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*we(i, 1, k)
                        end do
                    end do
                end do

                do k=1, nt
                    do i=1, imm1
                        do np1=3, odd_stride, 2
                            br(1, np1, k) = br(1, np1, k)+vb(i, np1)*vo(i, 1, k)
                        end do
                    end do
                end do
                !
                !  case m = 1 through nlat-1
                !
                if (mmax < 2) return

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*nlat-(m*(m + 1))/2
                    mp2 = mp1+1

                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do i=1, imm1
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*we(i, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*we(i, 2*mp1-2, k)
                                end do
                            end do
                        end do
                        if (mod(nlat, 2) /= 0) then
                            do k=1, nt
                                do np1=mp1, odd_stride, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+wb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)-wb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                                end do
                            end do
                        end if
                    end if

                    if (mp2 > even_stride) return

                    do k=1, nt
                        do i=1, imm1
                            do np1=mp2, even_stride, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-2, k) &
                                    +wb(i, np1+mb)*vo(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-1, k) &
                                    -wb(i, np1+mb)*vo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do

                    if (mod(nlat, 2) == 0) return

                    do k=1, nt
                        do np1=mp2, even_stride, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-1, k)
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
                            br(1, np1, k) = br(1, np1, k)+vb(i, np1)*vo(i, 1, k)
                        end do
                    end do
                end do
                !
                !  case m = 1 through nlat-1
                !
                if (mmax < 2) return

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*nlat-(m*(m + 1))/2
                    mp2 = mp1+1

                    if (mp1 > odd_stride) exit

                    do k=1, nt
                        do i=1, imm1
                            do np1=mp1, odd_stride, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-2, k) &
                                    +wb(i, np1+mb)*we(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-1, k) &
                                    -wb(i, np1+mb)*we(i, 2*mp1-2, k)
                            end do
                        end do
                    end do


                    if (mod(nlat, 2) == 0) return

                    do  k=1, nt
                        do np1=mp1, odd_stride, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+wb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)-wb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                        end do
                    end do
                end do
            case (8)
                !
                !  case ityp=8   v odd, w even, and both br and bi equal zero
                !
                ! case m=0
                !
                do k=1, nt
                    do i=1, imid
                        do np1=2, even_stride, 2
                            cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*we(i, 1, k)
                        end do
                    end do
                end do
                !
                !  case m = 1 through nlat-1
                !
                if (mmax < 2) return

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*nlat-(m*(m + 1))/2
                    mp2 = mp1+1

                    if (mp2 > even_stride) return

                    do k=1, nt
                        do i=1, imm1
                            do np1=mp2, even_stride, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-2, k) &
                                    +wb(i, np1+mb)*vo(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-1, k) &
                                    -wb(i, np1+mb)*vo(i, 2*mp1-2, k)
                            end do
                        end do
                    end do

                    if (mod(nlat, 2) == 0) return

                    do  k=1, nt
                        do  np1=mp2, even_stride, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                        end do
                    end do
                end do
        end select vector_symmetry_cases

    end subroutine vhags_lower_utility_routine

    pure function get_vhagsi_workspace_indices(nlat, imid, lmn) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: imid
        integer(ip), intent(in)  :: lmn
        integer(ip)              :: return_value(7)

        ! Set workspace pointers for indices
        associate (i => return_value)
            i(1) = 1
            i(2) = i(1)+imid*lmn
            i(3) = i(2)+imid*lmn
            i(4) = 1
            i(5) = i(4)+ nlat
            i(6) = i(5)+ nlat
            i(7) = i(6)+3*imid*nlat
        end associate

    end function get_vhagsi_workspace_indices

    ! Purpose:
    !
    ! Computes associated legendre functions.
    subroutine precompute_associated_legendre_functions(nlat, imid, vb, wb, dthet, dwts, dpbar, work)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: imid
        real(wp),    intent(out) :: vb(imid, *)
        real(wp),    intent(out) :: wb(imid, *)
        real(wp),    intent(out) :: dthet(nlat)
        real(wp),    intent(out) :: dwts(nlat)
        real(wp),    intent(out) :: dpbar(imid, nlat, 3)
        real(wp),    intent(out) :: work(*)

        ! Local variables
        integer(ip)         :: i, local_error_flag, id, ix, iy
        integer(ip)         :: m, mn, n, nm, np, nz
        real(wp)            :: abel, bbel, cbel, dcf
        type(SpherepackUtility) :: util

        !  Compute gaussian grid
        call compute_gaussian_latitudes_and_weights(nlat, dthet, dwts, local_error_flag)

        ! Set m=n=0 legendre polynomials for all theta(i)
        dpbar(:, 1, 1) = cos(PI/4)
        vb(:, 1) = ZERO
        wb(:, 1) = ZERO

        ! Main loop for remaining vb, and wb
        main_loop: do n=1, nlat - 1

            nm = mod(n - 2, 3) + 1
            nz = mod(n - 1, 3) + 1
            np = mod(n, 3) + 1

            !  Compute dpbar for m=0
            call util%compute_fourier_coefficients(0, n, work)
            mn = get_index(0, n, nlat)
            do i=1, imid
                call util%compute_legendre_polys_from_fourier_coeff(0, n, dthet(i), work, dpbar(i, 1, np))
            end do

            !  Compute dpbar for m=1
            call util%compute_fourier_coefficients(1, n, work)

            mn = get_index(1, n, nlat)

            do i=1, imid
                call util%compute_legendre_polys_from_fourier_coeff(1, n, dthet(i), work, dpbar(i, 2, np))
            end do

            ! Compute and store dpbar for m=2, n
            if (2 <= n) then
                do m=2, n
                    abel = sqrt(real((2*n + 1)*(m + n-2)*(m + n - 3))/ &
                        real((2*n - 3)*(m + n - 1)*(m + n)))
                    bbel = sqrt(real((2*n + 1)*(n - m - 1)*(n - m))/ &
                        real((2*n - 3)*(m + n - 1)*(m + n)))
                    cbel = sqrt(real((n - m + 1)*(n - m + 2))/ &
                        real((m + n - 1)*(m + n)))
                    id = get_index(m, n, nlat)

                    if (n - 1 <= m) then
                        dpbar(1:imid, m + 1, np) = &
                            abel*dpbar(1:imid, m - 1, nm)-cbel*dpbar(1:imid, m - 1, np)
                    else
                        dpbar(1:imid, m + 1, np) = &
                            abel*dpbar(1:imid, m - 1, nm)+bbel*dpbar(1:imid, m + 1, nm) &
                            -cbel*dpbar(1:imid, m - 1, np)
                    end if
                end do
            end if

            ! Compute the derivative of the functions
            ix = get_index(0, n, nlat)
            iy = get_index(n, n, nlat)
            vb(1:imid, ix) = -dpbar(1:imid, 2, np)*dwts(1:imid)
            vb(1:imid, iy) = dpbar(1:imid, n, np)/sqrt(real(2*(n + 1), kind=wp))*dwts(1:imid)

            select case (n)
                case (1)
                    ! Compute the vector harmonic w(theta) = m*pbar/cos(theta)
                    !
                    ! Set wb=0 for m=0
                    ix = get_index(0, n, nlat)
                    wb(1:imid, ix) = ZERO
                case default
                    dcf = sqrt(real(4*n*(n + 1), kind=wp))
                    do m=1, n - 1
                        ix = get_index(m, n, nlat)
                        abel = sqrt(real((n + m)*(n - m + 1), kind=wp))/dcf
                        bbel = sqrt(real((n - m)*(n + m + 1), kind=wp))/dcf
                        vb(1:imid, ix) = &
                            (abel * dpbar(1:imid, m, np) - bbel * dpbar(1:imid, m + 2, np))&
                            * dwts(1:imid)
                    end do
            end select

            ! Compute wb for m=1, n
            dcf = sqrt(real(2*n + 1, kind=wp)/real(4*n*(n + 1)*(n + n - 1), kind=wp))
            do m=1, n
                ix = get_index(m, n, nlat)
                abel = dcf * sqrt(real((n + m)*(n + m - 1), kind=wp))
                bbel = dcf * sqrt(real((n - m)*(n - m - 1), kind=wp))
                if (n - 1 <= m) then
                    wb(1:imid, ix) = abel * dpbar(1:imid, m, nz) * dwts(1:imid)
                else
                    wb(1:imid, ix) = &
                        (abel * dpbar(1:imid, m, nz) + bbel * dpbar(1:imid, m + 2, nz))&
                        * dwts(1:imid)
                end if
            end do
        end do main_loop

    end subroutine precompute_associated_legendre_functions

    pure function get_index(m, n, nlat) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: m
        integer(ip), intent(in) :: n
        integer(ip), intent(in) :: nlat
        integer(ip)             :: return_value

        return_value = m*nlat-(m*(m + 1))/2 + n + 1

    end function get_index

end submodule vector_analysis_gaussian_grid_saved
