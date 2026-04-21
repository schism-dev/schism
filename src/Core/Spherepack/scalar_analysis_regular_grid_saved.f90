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
submodule(scalar_analysis_routines) scalar_analysis_regular_grid_saved

! Parameter confined to the module
integer(ip), parameter :: NUMBER_OF_WORKSPACE_INDICES = 3

contains

    ! Purpose:
    !
    !     subroutine shaes(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, wshaes, ierror)
    !
    !     subroutine shaes performs the spherical harmonic analysis
    !     on the array g and stores the result in the arrays a and b.
    !     the analysis is performed on an equally spaced grid.  the
    !     associated legendre functions are stored rather than recomputed
    !     as they are in subroutine shaec.  the analysis is described
    !     below at output parameters a, b.
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
    !            than or equal to 4. the efficiency of the computation is
    !            improved when nlon is a product of small prime numbers.
    !
    !     isym   = 0  no symmetries exist about the equator. the analysis
    !                 is performed on the entire sphere.  i.e. on the
    !                 array g(i, j) for i=1, ..., nlat and j=1, ..., nlon.
    !                 (see description of g below)
    !
    !            = 1  g is antisymmetric about the equator. the analysis
    !                 is performed on the northern hemisphere only.  i.e.
    !                 if nlat is odd the analysis is performed on the
    !                 array g(i, j) for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon.
    !                 if nlat is even the analysis is performed on the
    !                 array g(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !
    !            = 2  g is symmetric about the equator. the analysis is
    !                 performed on the northern hemisphere only.  i.e.
    !                 if nlat is odd the analysis is performed on the
    !                 array g(i, j) for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon.
    !                 if nlat is even the analysis is performed on the
    !                 array g(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !     nt     the number of analyses.  in the program that calls shaes, 
    !            the arrays g, a and b can be three dimensional in which
    !            case multiple analyses will be performed.  the third
    !            index is the analysis index which assumes the values
    !            k=1, ..., nt.  for a single analysis set nt=1. the
    !            discription of the remaining parameters is simplified
    !            by assuming that nt=1 or that the arrays g, a and b
    !            have only two dimensions.
    !
    !     g      a two or three dimensional array (see input parameter
    !            nt) that contains the discrete function to be analyzed.
    !            g(i, j) contains the value of the function at the colatitude
    !            point theta(i) = (i-1)*pi/(nlat-1) and longitude point
    !            phi(j) = (j-1)*2*pi/nlon. the index ranges are defined
    !            above at the input parameter isym.
    !
    !
    !     idg    the first dimension of the array g as it appears in the
    !            program that calls shaes.  if isym equals zero then idg
    !            must be at least nlat.  if isym is nonzero then idg
    !            must be at least nlat/2 if nlat is even or at least
    !            (nlat + 1)/2 if nlat is odd.
    !
    !     jdg    the second dimension of the array g as it appears in the
    !            program that calls shaes.  jdg must be at least nlon.
    !
    !     mdab   the first dimension of the arrays a and b as it appears
    !            in the program that calls shaes. mdab must be at least
    !            min(nlat, (nlon+2)/2) if nlon is even or at least
    !            min(nlat, (nlon + 1)/2) if nlon is odd.
    !
    !     ndab   the second dimension of the arrays a and b as it appears
    !            in the program that calls shaes. ndab must be at least nlat
    !
    !     wshaes an array which must be initialized by subroutine shaesi.
    !            once initialized, wshaes can be used repeatedly by shaes
    !            as long as nlon and nlat remain unchanged.  wshaes must
    !            not be altered between calls of shaes.
    !
    !     lshaes the dimension of the array wshaes as it appears in the
    !            program that calls shaes. define
    !
    !               l1 = min(nlat, (nlon+2)/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lshaes must be at least
    !
    !               (l1*l2*(2*nlat-l1+1))/2+nlon+15
    !
    !     output parameters
    !
    !     a, b    both a, b are two or three dimensional arrays (see input
    !            parameter nt) that contain the spherical harmonic
    !            coefficients in the representation of g(i, j) given in the
    !            discription of subroutine shses. for isym=0, a(m, n) and
    !            b(m, n) are given by the equations listed below. symmetric
    !            versions are used when isym is greater than zero.
    !
    !
    !
    !     definitions
    !
    !     1. the normalized associated legendre functions
    !
    !     pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
    !                       *sin(theta)**m/(2**n*factorial(n)) times the
    !                       (n+m)th derivative of (x**2-1)**n with respect
    !                       to x=cos(theta)
    !
    !     2. the normalized z functions for m even
    !
    !     zbar(m, n, theta) = 2/(nlat-1) times the sum from k=0 to k=nlat-1 of
    !                       the integral from tau = 0 to tau = pi of
    !                       cos(k*theta)*cos(k*tau)*pbar(m, n, tau)*sin(tau)
    !                       (first and last terms in this sum are divided
    !                       by 2)
    !
    !     3. the normalized z functions for m odd
    !
    !     zbar(m, n, theta) = 2/(nlat-1) times the sum from k=0 to k=nlat-1 of
    !                       of the integral from tau = 0 to tau = pi of
    !                       sin(k*theta)*sin(k*tau)*pbar(m, n, tau)*sin(tau)
    !
    !     4. the fourier transform of g(i, j).
    !
    !     c(m, i)          = 2/nlon times the sum from j=1 to j=nlon
    !                       of g(i, j)*cos((m-1)*(j-1)*2*pi/nlon)
    !                       (the first and last terms in this sum
    !                       are divided by 2)
    !
    !     s(m, i)          = 2/nlon times the sum from j=2 to j=nlon
    !                       of g(i, j)*sin((m-1)*(j-1)*2*pi/nlon)
    !
    !     5. the maximum (plus one) longitudinal wave number
    !
    !            mmax = min(nlat, (nlon+2)/2) if nlon is even or
    !            mmax = min(nlat, (nlon + 1)/2) if nlon is odd.
    !
    !     then for m=0, ..., mmax-1 and  n=m, ..., nlat-1  the arrays a, b are
    !     given by
    !
    !     a(m+1, n+1)      = the sum from i=1 to i=nlat of
    !                       c(m+1, i)*zbar(m, n, theta(i))
    !                       (first and last terms in this sum are
    !                       divided by 2)
    !
    !     b(m+1, n+1)      = the sum from i=1 to i=nlat of
    !                       s(m+1, i)*zbar(m, n, theta(i))
    !
    !
    !     ierror = 0  no errors
    !            = 1  error in the specification of nlat
    !            = 2  error in the specification of nlon
    !            = 3  error in the specification of isym
    !            = 4  error in the specification of nt
    !            = 5  error in the specification of idg
    !            = 6  error in the specification of jdg
    !            = 7  error in the specification of mdab
    !            = 8  error in the specification of ndab
    !            = 9  error in the specification of lshaes
    !
    module subroutine shaes(nlat, nlon, isym, nt, g, idg, jdg, a, b, &
        mdab, ndab, wshaes, ierror)

        ! Dummy arguments
        integer(ip), intent(in)   :: nlat
        integer(ip), intent(in)   :: nlon
        integer(ip), intent(in)   :: isym
        integer(ip), intent(in)   :: nt
        real(wp),    intent(in)   :: g(idg, jdg, nt)
        integer(ip), intent(in)   :: idg
        integer(ip), intent(in)   :: jdg
        real(wp),    intent(out)  :: a(mdab, ndab, nt)
        real(wp),    intent(out)  :: b(mdab, ndab, nt)
        integer(ip), intent(in)   :: mdab
        integer(ip), intent(in)   :: ndab
        real(wp),    intent(in)   :: wshaes(:)
        integer(ip), intent(out)  :: ierror

        ! Local variables
        integer(ip) :: mmax, imid, idz, lzimn, ls, ifft, ist, lwork, nln
        integer(ip) :: required_wavetable_size
        type(SpherepackUtility) :: util

        ! Check input arguments
        required_wavetable_size = util%get_lshaes(nlat, nlon)

        call util%check_scalar_transform_inputs(isym, idg, jdg, &
            mdab, ndab, nlat, nlon, nt, required_wavetable_size, &
            wshaes, ierror)

        ! Check error flag
        if (ierror /= 0) return

        mmax = min(nlat, nlon/2+1)
        imid = (nlat + 1)/2
        idz = (mmax*(2*nlat-mmax+1))/2
        lzimn = idz*imid
        ifft = lzimn + 1

        !  Set calling argument for analysis
        select case (isym)
            case (0)
                ls = nlat
                ist = imid
            case default
                ls = imid
                ist = 0
        end select

        nln = nt*ls*nlon

        ! Set required workspace size
        lwork = nln+ls*nlon

        block
            integer(ip) :: iw1
            real(wp)    :: work(lwork)

            ! Set workspace pointer indices
            iw1 = ist + 1

            !  Perform analysis
            call shaes_lower_utility_routine(nlat, isym, nt, g, idg, jdg, a, b, &
                mdab, ndab, wshaes, idz, ls, nlon, work, work(iw1:), wshaes(ifft:))
        end block

    end subroutine shaes

    ! Purpose:
    !
    !     subroutine shaesi(nlat, nlon, wshaes, lshaes, work, lwork, dwork, ldwork, ierror)
    !
    !     subroutine shaesi initializes the array wshaes which can then
    !     be used repeatedly by subroutine shaes
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
    !            than or equal to 4. the efficiency of the computation is
    !            improved when nlon is a product of small prime numbers.
    !
    !     lshaes the dimension of the array wshaes as it appears in the
    !            program that calls shaesi. define
    !
    !               l1 = min(nlat, (nlon+2)/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lshaes must be at least
    !
    !               (l1*l2*(2*nlat-l1+1))/2+nlon+15
    !
    !     work   a real   work array that does not have to be saved.
    !
    !     lwork  the dimension of the array work as it appears in the
    !            program that calls shaesi.  define
    !
    !               l1 = min(nlat, (nlon+2)/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lwork must be at least
    !
    !               5*nlat*l2+3*((l1-2)*(2*nlat-l1-1))/2
    !
    !
    !     dwork  a real work array that does not have to be saved.
    !
    !     ldwork the dimension of the array dwork as it appears in the
    !            program that calls shaesi.  ldwork must be at least nlat+1
    !
    !
    !     output parameters
    !
    !     wshaes an array which is initialized for use by subroutine shaes.
    !            once initialized, wshaes can be used repeatedly by shaes
    !            as long as nlon and nlat remain unchanged.  wshaes must
    !            not be altered between calls of shaes.
    !
    !     ierror = 0  no errors
    !            = 1  error in the specification of nlat
    !            = 2  error in the specification of nlon
    !            = 3  error in the specification of lshaes
    !            = 4  error in the specification of lwork
    !            = 5  error in the specification of ldwork
    !
    ! Remarks:
    !
    ! size(wshaes) = (l*(l+1)*imid)/2+nlon+15
    ! size(work) = 5*l*imid + 3*((l-3)*l+2)/2
    !
    module subroutine shaesi(nlat, nlon, wshaes, ierror)

        ! Dummy arguments
        integer(ip), intent(in)    :: nlat
        integer(ip), intent(in)    :: nlon
        real(wp),    intent(out)   :: wshaes(:)
        integer(ip), intent(out)   :: ierror

        ! Local variables
        integer(ip) :: mmax, imid, labc, lzimn
        integer(ip) :: workspace_indices(NUMBER_OF_WORKSPACE_INDICES)
        integer(ip) :: lwork, required_wavetable_size
        type(SpherepackUtility) :: util

        mmax = min(nlat, nlon/2+1)
        imid = (nlat + 1)/2
        labc = 3*((mmax-2)*(2*nlat-mmax-1))/2
        lzimn = (imid*mmax*(2*nlat-mmax+1))/2
        required_wavetable_size = lzimn+nlon+15

        ! Check calling arguments
        if (nlat < 3) then
            ierror = 1
        else if (nlon < 4) then
            ierror = 2
        else if (size(wshaes) < required_wavetable_size) then
            ierror = 3
        else
            ierror = 0
        end if

        ! Address error flag
        if (ierror /= 0) return

        ! Set required workspace sizes
        lwork = 5*nlat*imid + labc

        !  Set workspace indices
        workspace_indices = get_workspace_indices(nlat, nlon, mmax, imid, lzimn)

        !  Compute wavetable
        block
            real(wp) :: work(lwork)
            associate (&
                idz => workspace_indices(1), &
                iw1 => workspace_indices(2), &
                iw2 => workspace_indices(3) &
                )
                call util%initialize_scalar_analysis_regular_grid_saved( &
                    nlat, nlon, imid, wshaes, idz, work, work(iw1:))

                call util%hfft%initialize(nlon, wshaes(iw2:))
            end associate
        end block

    end subroutine shaesi

    subroutine shaes_lower_utility_routine(nlat, isym, nt, g, idgs, jdgs, &
        a, b, mdab, ndab, z, idz, idg, jdg, g_even, g_odd, whrfft)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: isym
        integer(ip), intent(in)  :: nt
        real(wp),    intent(in)  :: g(idgs, jdgs, nt)
        integer(ip), intent(in)  :: idgs
        integer(ip), intent(in)  :: jdgs
        real(wp),    intent(out) :: a(mdab, ndab, nt)
        real(wp),    intent(out) :: b(mdab, ndab, nt)
        integer(ip), intent(in)  :: mdab
        integer(ip), intent(in)  :: ndab
        real(wp),    intent(in)  :: z(idz, *)
        integer(ip), intent(in)  :: idz
        integer(ip), intent(in)  :: idg
        integer(ip), intent(in)  :: jdg
        real(wp),    intent(out) :: g_even(idg, jdg, nt)
        real(wp),    intent(out) :: g_odd(idg, jdg, nt)
        real(wp),    intent(in)  :: whrfft(:)

        ! Local variables
        integer(ip)             :: i, k, m, mb, ls, mp1, np1, mp2
        integer(ip)             :: order_stride, degree_stride
        integer(ip)             :: imm1, imid, m_trunc, nlon
        real(wp)                :: fsn, tsn
        type(SpherepackUtility) :: util

        ls = idg
        nlon = jdg
        m_trunc = min(nlat, nlon/2+1)

        if (2*m_trunc-1 > nlon) then
            order_stride = m_trunc-1
        else
            order_stride = m_trunc
        end if

        tsn = TWO/nlon
        fsn = FOUR/nlon
        imid = (nlat + 1)/2

        if (even(nlat)) then
            imm1 = imid
        else
            imm1 = imid-1
        end if

        select case (isym)
            case(0)
                do i=1, imm1
                    g_even(i, 1:nlon, :) = tsn * (g(i, 1:nlon, :) + g((nlat + 1)-i, 1:nlon, :))
                    g_odd(i, 1:nlon, :) = tsn * (g(i, 1:nlon, :) - g((nlat + 1)-i, 1:nlon, :))
                end do
            case default
                g_even(1:imm1, 1:nlon, :) = fsn * g(1:imm1, 1:nlon, :)
        end select

        if ((isym /= 1) .and. odd(nlat)) then
            g_even(imid, 1:nlon, :) = tsn * g(imid, 1:nlon, :)
        end if

        !  Fast Fourier Transform
        fft_loop: do k=1, nt

            call util%hfft%forward(ls, nlon, g_even(:, :, k), ls, whrfft)

            if (odd(nlon)) exit fft_loop

            g_even(1:ls, nlon, k) = HALF * g_even(1:ls, nlon, k)

        end do fft_loop

        ! Preset coefficients to zero
        a = ZERO
        b = ZERO

        if (isym /= 1) then

            do k=1, nt
                do i=1, imid
                    do np1=1, nlat, 2
                        a(1, np1, k) = a(1, np1, k)+z(np1, i)*g_even(i, 1, k)
                    end do
                end do
            end do

            if (even(nlat)) then
                degree_stride = nlat-1
            else
                degree_stride = nlat
            end if

            do mp1=2, order_stride
                m = mp1-1
                mb = m*(nlat-1)-(m*(m-1))/2
                do k=1, nt
                    do i=1, imid
                        do np1=mp1, degree_stride, 2
                            a(mp1, np1, k) = a(mp1, np1, k)+z(np1+mb, i)*g_even(i, 2*mp1-2, k)
                            b(mp1, np1, k) = b(mp1, np1, k)+z(np1+mb, i)*g_even(i, 2*mp1-1, k)
                        end do
                    end do
                end do
            end do

            if (order_stride /= m_trunc .and. m_trunc <= degree_stride) then
                mb = order_stride*(nlat-1)-(order_stride*(order_stride-1))/2
                do k=1, nt
                    do i=1, imid
                        do np1=m_trunc, degree_stride, 2
                            a(m_trunc, np1, k) = a(m_trunc, np1, k)+z(np1+mb, i)*g_even(i, 2*m_trunc-2, k)
                        end do
                    end do
                end do
            end if
            if (isym == 2) return
        end if

        do k=1, nt
            do i=1, imm1
                do np1=2, nlat, 2
                    a(1, np1, k) = a(1, np1, k)+z(np1, i)*g_odd(i, 1, k)
                end do
            end do
        end do

        if (even(nlat)) then
            degree_stride = nlat
        else
            degree_stride = nlat-1
        end if

        do mp1=2, order_stride
            m = mp1-1
            mp2 = mp1+1
            mb = m*(nlat-1)-(m*(m-1))/2
            do k=1, nt
                do i=1, imm1
                    do np1=mp2, degree_stride, 2
                        a(mp1, np1, k) = a(mp1, np1, k)+z(np1+mb, i)*g_odd(i, 2*mp1-2, k)
                        b(mp1, np1, k) = b(mp1, np1, k)+z(np1+mb, i)*g_odd(i, 2*mp1-1, k)
                    end do
                end do
            end do
        end do

        mp2 = m_trunc+1
        if (order_stride /= m_trunc .and. mp2 <= degree_stride) then
            mb = order_stride*(nlat-1)-(order_stride*(order_stride-1))/2
            do k=1, nt
                do i=1, imm1
                    do np1=mp2, degree_stride, 2
                        a(m_trunc, np1, k) = a(m_trunc, np1, k)+z(np1+mb, i)*g_odd(i, 2*m_trunc-2, k)
                    end do
                end do
            end do
        end if

    end subroutine shaes_lower_utility_routine

    pure function get_workspace_indices(nlat, nlon, mmax, imid, lzimn) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip), intent(in) :: mmax
        integer(ip), intent(in) :: imid
        integer(ip), intent(in) :: lzimn
        integer(ip)             :: return_value(NUMBER_OF_WORKSPACE_INDICES)

        associate (i => return_value)
            i(1) = (mmax*(2*nlat-mmax+1))/2
            i(2) = 3*nlat*imid+1
            i(3) = lzimn + 1
        end associate

    end function get_workspace_indices

end submodule scalar_analysis_regular_grid_saved
