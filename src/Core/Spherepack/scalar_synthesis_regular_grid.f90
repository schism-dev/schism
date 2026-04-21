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
submodule(scalar_synthesis_routines) scalar_synthesis_regular_grid

contains

    !
    !     subroutine shsec(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
    !                      wshsec, ierror)
    !
    !     subroutine shsec performs the spherical harmonic synthesis
    !     on the arrays a and b and stores the result in the array g.
    !     the synthesis is performed on an equally spaced grid.  the
    !     associated legendre functions are recomputed rather than stored
    !     as they are in subroutine shses.  the synthesis is described
    !     below at output parameter g.
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
    !     isym   = 0  no symmetries exist about the equator. the synthesis
    !                 is performed on the entire sphere.  i.e. on the
    !                 array g(i, j) for i=1, ..., nlat and j=1, ..., nlon.
    !                 (see description of g below)
    !
    !            = 1  g is antisymmetric about the equator. the synthesis
    !                 is performed on the northern hemisphere only.  i.e.
    !                 if nlat is odd the synthesis is performed on the
    !                 array g(i, j) for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon.
    !                 if nlat is even the synthesis is performed on the
    !                 array g(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !
    !            = 2  g is symmetric about the equator. the synthesis is
    !                 performed on the northern hemisphere only.  i.e.
    !                 if nlat is odd the synthesis is performed on the
    !                 array g(i, j) for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon.
    !                 if nlat is even the synthesis is performed on the
    !                 array g(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !     nt     the number of syntheses.  in the program that calls shsec,
    !            the arrays g, a and b can be three dimensional in which
    !            case multiple syntheses will be performed.  the third
    !            index is the synthesis index which assumes the values
    !            k=1, ..., nt.  for a single synthesis set nt=1. the
    !            discription of the remaining parameters is simplified
    !            by assuming that nt=1 or that the arrays g, a and b
    !            have only two dimensions.
    !
    !     idg    the first dimension of the array g as it appears in the
    !            program that calls shsec.  if isym equals zero then idg
    !            must be at least nlat.  if isym is nonzero then idg
    !            must be at least nlat/2 if nlat is even or at least
    !            (nlat + 1)/2 if nlat is odd.
    !
    !     jdg    the second dimension of the array g as it appears in the
    !            program that calls shsec.  jdg must be at least nlon.
    !
    !     a, b    two or three dimensional arrays (see the input parameter
    !            nt) that contain the coefficients in the spherical harmonic
    !            expansion of g(i, j) given below at the definition of the
    !            output parameter g.  a(m, n) and b(m, n) are defined for
    !            indices m=1, ..., mmax and n=m, ..., nlat where mmax is the
    !            maximum (plus one) longitudinal wave number given by
    !            mmax = min(nlat, (nlon+2)/2) if nlon is even or
    !            mmax = min(nlat, (nlon + 1)/2) if nlon is odd.
    !
    !     mdab   the first dimension of the arrays a and b as it appears
    !            in the program that calls shsec. mdab must be at least
    !            min(nlat, (nlon+2)/2) if nlon is even or at least
    !            min(nlat, (nlon + 1)/2) if nlon is odd.
    !
    !     ndab   the second dimension of the arrays a and b as it appears
    !            in the program that calls shsec. ndab must be at least nlat
    !
    !     wshsec an array which must be initialized by subroutine shseci.
    !            once initialized, wshsec can be used repeatedly by shsec
    !            as long as nlon and nlat remain unchanged.  wshsec must
    !            not be altered between calls of shsec.
    !
    !     lshsec the dimension of the array wshsec as it appears in the
    !            program that calls shsec. define
    !
    !               l1 = min(nlat, (nlon+2)/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lshsec must be at least
    !
    !            2*nlat*l2+3*((l1-2)*(2*nlat-l1-1))/2+nlon+15
    !
    !     output parameters
    !
    !     g      a two or three dimensional array (see input parameter
    !            nt) that contains the spherical harmonic synthesis of
    !            the arrays a and b at the colatitude point theta(i) =
    !            (i-1)*pi/(nlat-1) and longitude point phi(j) =
    !            (j-1)*2*pi/nlon. the index ranges are defined above at
    !            at the input parameter isym.  for isym=0, g(i, j) is
    !            given by the the equations listed below.  symmetric
    !            versions are used when isym is greater than zero.
    !
    !     the normalized associated legendre functions are given by
    !
    !     pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
    !                       *sin(theta)**m/(2**n*factorial(n)) times the
    !                       (n+m)th derivative of (x**2-1)**n with respect
    !                       to x=cos(theta)
    !
    !     define the maximum (plus one) longitudinal wave number
    !     as   mmax = min(nlat, (nlon+2)/2) if nlon is even or
    !          mmax = min(nlat, (nlon + 1)/2) if nlon is odd.
    !
    !     then g(i, j) = the sum from n=0 to n=nlat-1 of
    !
    !                   HALF*pbar(0, n, theta(i))*a(1, n+1)
    !
    !              plus the sum from m=1 to m=mmax-1 of
    !
    !                   the sum from n=m to n=nlat-1 of
    !
    !              pbar(m, n, theta(i))*(a(m+1, n+1)*cos(m*phi(j))
    !                                    -b(m+1, n+1)*sin(m*phi(j)))
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
    !            = 9  error in the specification of lshsec
    !
    module subroutine shsec(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, wshsec, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: isym
        integer(ip), intent(in)  :: nt
        real(wp),    intent(out) :: g(idg, jdg, nt)
        integer(ip), intent(in)  :: idg
        integer(ip), intent(in)  :: jdg
        real(wp),    intent(in)  :: a(mdab, ndab, nt)
        real(wp),    intent(in)  :: b(mdab, ndab, nt)
        integer(ip), intent(in)  :: mdab
        integer(ip), intent(in)  :: ndab
        real(wp),    intent(in)  :: wshsec(:)
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip) :: imid, ist, labc, ls
        integer(ip) :: lzz1, mmax, nln, lwork
        integer(ip)             :: required_wavetable_size
        type(SpherepackUtility) :: util

        ! Check input arguments
        required_wavetable_size = util%get_lshsec(nlat, nlon)

        call util%check_scalar_transform_inputs(isym, idg, jdg, &
            mdab, ndab, nlat, nlon, nt, required_wavetable_size, &
            wshsec, ierror)

        ! Check error flag
        if (ierror /= 0) return

        mmax = min(nlat, nlon/2+1)
        imid = (nlat + 1)/2
        lzz1 = 2*nlat*imid
        labc = 3*((mmax-2)*(2*nlat-mmax-1))/2

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
        lwork = nln+max(ls*nlon, 3*nlat*imid)

        block
            real(wp)    :: work(lwork)
            integer(ip) :: jw1, jw2, iw1

            jw1 = ist + 1
            jw2 = nln + 1
            iw1 = lzz1+labc+1

            call shsec_lower_utility_routine(nlat, isym, nt, g, idg, jdg, a, b, &
                mdab, ndab, imid, ls, nlon, work, work(jw1:), work(jw2:), work(jw2:), wshsec, wshsec(iw1:))
        end block

    end subroutine shsec

    !     subroutine shseci(nlat, nlon, wshsec, ierror)
    !
    !     subroutine shseci initializes the array wshsec which can then
    !     be used repeatedly by subroutine shsec.
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
    !     lshsec the dimension of the array wshsec as it appears in the
    !            program that calls shseci. the array wshsec is an output
    !            parameter which is described below. define
    !
    !               l1 = min(nlat, (nlon+2)/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lshsec must be at least
    !
    !            2*nlat*l2+3*((l1-2)*(2*nlat-l1-1))/2+nlon+15
    !
    !     output parameters
    !
    !     wshsec an array which is initialized for use by subroutine shsec.
    !            once initialized, wshsec can be used repeatedly by shsec
    !            as long as nlon and nlat remain unchanged.  wshsec must
    !            not be altered between calls of shsec.
    !
    !     ierror = 0  no errors
    !            = 1  error in the specification of nlat
    !            = 2  error in the specification of nlon
    !            = 3  error in the specification of lshsec
    !
    module subroutine shseci(nlat, nlon, wshsec, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(out) :: wshsec(:)
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip) :: imid, iw1
        integer(ip) :: labc, lzz1, mmax, ldwork
        type(SpherepackUtility) :: util

        associate (lshsec => size(wshsec))

            imid = (nlat + 1)/2
            mmax = min(nlat, nlon/2+1)
            lzz1 = 2*nlat*imid
            labc = 3*((mmax-2)*(2*nlat-mmax-1))/2

            ! Check calling arguments
            if (nlat < 3) then
                ierror = 1
            else if (nlon < 4) then
                ierror = 2
            else if (lshsec < lzz1+labc+nlon+15) then
                ierror = 3
            else
                ierror = 0
            end if

            ! Check error flag
            if (ierror /= 0) return

            ! Set required workspace size
            ldwork = nlat + 1

            block
                real(wp) :: dwork(ldwork)
                call util%initialize_scalar_synthesis_regular_grid(nlat, nlon, wshsec, dwork)
            end block

            ! set workspace index
            iw1 = lzz1+labc+1

            call util%hfft%initialize(nlon, wshsec(iw1:))
        end associate

    end subroutine shseci

    subroutine shsec_lower_utility_routine(nlat, isym, nt, g, idgs, jdgs, a, b, mdab, ndab, imid, &
        idg, jdg, ge, go, work, pb, walin, whrfft)

        real(wp) :: a
        real(wp) :: b
        real(wp) :: g
        real(wp) :: ge
        real(wp) :: go
        integer(ip) :: i
        integer(ip) :: i3
        integer(ip) :: idg
        integer(ip) :: idgs
        integer(ip) :: imid
        integer(ip) :: imm1
        integer(ip) :: isym
        integer(ip) :: j
        integer(ip) :: jdg
        integer(ip) :: jdgs
        integer(ip) :: k
        integer(ip) :: ls
        integer(ip) :: m
        integer(ip) :: mdab
        integer(ip) :: mdo
        integer(ip) :: mmax
        integer(ip) :: modl
        integer(ip) :: mp1
        integer(ip) :: mp2
        integer(ip) :: ndab
        integer(ip) :: ndo
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: nlp1
        integer(ip) :: np1
        integer(ip) :: nt
        real(wp) :: pb
        real(wp) :: walin
        real(wp) :: whrfft
        real(wp) :: work


        ! Purpose:
        !
        ! Requirements
        !
        ! whrfft must have at least nlon+15 locations
        ! walin must have 3*l*imid + 3*((l-3)*l+2)/2 locations
        ! zb must have 3*l*imid locations
        !
        dimension g(idgs, jdgs, nt), a(mdab, ndab, nt), b(mdab, ndab, nt), &
            ge(idg, jdg, nt), go(idg, jdg, nt), pb(imid, nlat, 3), walin(*), &
            whrfft(*), work(*)

        type(SpherepackUtility) :: util

        ls = idg
        nlon = jdg
        mmax = min(nlat, nlon/2+1)

        if (2*mmax-1 > nlon) then
            mdo = mmax-1
        else
            mdo = mmax
        end if

        nlp1 = nlat+1
        modl = mod(nlat, 2)

        if (modl /= 0) then
            imm1 = imid-1
        else
            imm1 = imid
        end if

        do k=1, nt
            ge(1: ls, 1: nlon, k)=0.0_wp
        end do

        if_block: block

            if (isym /= 1) then
                call util%compute_legendre_polys_regular_grid(2, nlat, nlon, 0, pb, i3, walin)
                do k=1, nt
                    do np1=1, nlat, 2
                        do i=1, imid
                            ge(i, 1, k)=ge(i, 1, k)+a(1, np1, k)*pb(i, np1, i3)
                        end do
                    end do
                end do

                if (mod(nlat, 2) == 0) then
                    ndo = nlat-1
                else
                    ndo = nlat
                end if

                do mp1=2, mdo
                    m = mp1-1
                    call util%compute_legendre_polys_regular_grid(2, nlat, nlon, m, pb, i3, walin)
                    do np1=mp1, ndo, 2
                        do k=1, nt
                            do i=1, imid
                                ge(i, 2*mp1-2, k) = ge(i, 2*mp1-2, k)+a(mp1, np1, k)*pb(i, np1, i3)
                                ge(i, 2*mp1-1, k) = ge(i, 2*mp1-1, k)+b(mp1, np1, k)*pb(i, np1, i3)
                            end do
                        end do
                    end do
                end do

                if (.not.(mdo == mmax .or. mmax > ndo)) then
                    call util%compute_legendre_polys_regular_grid(2, nlat, nlon, mdo, pb, i3, walin)
                    do np1=mmax, ndo, 2
                        do k=1, nt
                            do i=1, imid
                                ge(i, 2*mmax-2, k) = ge(i, 2*mmax-2, k)+a(mmax, np1, k)*pb(i, np1, i3)
                            end do
                        end do
                    end do
                end if
                if (isym == 2) exit if_block
            end if

            call util%compute_legendre_polys_regular_grid(1, nlat, nlon, 0, pb, i3, walin)

            do k=1, nt
                do np1=2, nlat, 2
                    do i=1, imm1
                        go(i, 1, k)=go(i, 1, k)+a(1, np1, k)*pb(i, np1, i3)
                    end do
                end do
            end do


            if (mod(nlat, 2) /= 0) then
                ndo = nlat-1
            else
                ndo = nlat
            end if

            do mp1=2, mdo
                mp2 = mp1+1
                m = mp1-1
                call util%compute_legendre_polys_regular_grid(1, nlat, nlon, m, pb, i3, walin)
                do np1=mp2, ndo, 2
                    do k=1, nt
                        do i=1, imm1
                            go(i, 2*mp1-2, k) = go(i, 2*mp1-2, k)+a(mp1, np1, k)*pb(i, np1, i3)
                            go(i, 2*mp1-1, k) = go(i, 2*mp1-1, k)+b(mp1, np1, k)*pb(i, np1, i3)
                        end do
                    end do
                end do
            end do

            mp2 = mmax+1

            if (mdo == mmax .or. mp2 > ndo) exit if_block

            call util%compute_legendre_polys_regular_grid(1, nlat, nlon, mdo, pb, i3, walin)

            do np1=mp2, ndo, 2
                do k=1, nt
                    do i=1, imm1
                        go(i, 2*mmax-2, k) = go(i, 2*mmax-2, k)+a(mmax, np1, k)*pb(i, np1, i3)
                    end do
                end do
            end do

        end block if_block

        do k=1, nt

            if (mod(nlon, 2) == 0) ge(1: ls, nlon, k) = TWO*ge(1: ls, nlon, k)

            call util%hfft%backward(ls, nlon, ge(:, :, k), ls, whrfft)
        end do

        select case (isym)
            case (0)
                do k=1, nt
                    do j=1, nlon
                        do i=1, imm1
                            g(i, j, k) = HALF*(ge(i, j, k)+go(i, j, k))
                            g(nlp1-i, j, k) = HALF*(ge(i, j, k)-go(i, j, k))
                        end do
                        if (modl /= 0) g(imid, j, k) = HALF*ge(imid, j, k)
                    end do
                end do
            case default
                do k=1, nt
                    g(1: imid, 1: nlon, k) = HALF*ge(1: imid, 1: nlon, k)
                end do
        end select

    end subroutine shsec_lower_utility_routine

end submodule scalar_synthesis_regular_grid
