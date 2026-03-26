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
!     this file contains code and documentation for subroutines
!     shpgi and shpg.
!
submodule(scalar_projection_routines) scalar_projection_gaussian_grid

contains
    !
    ! Purpose:
    !
    !     shpg computes the harmonic projection, which is
    !     equivalent to a harmonic analysis (forward) followed
    !     by a harmonic synthesis (backward transform).
    !     shpg uses the n**2 projection or complement when appropriate
    !     as well as  odd/even factorization and zero truncation on an
    !     on a Gaussian distributed grid as defined in the JCP paper
    !     "Generalized discrete spherical harmonic transforms"
    !     by Paul N. Swarztrauber and William F. Spotz
    !     J. Comp. Phys., 159(2000) pp. 213-230.
    !
    !     subroutine shpg(nlat, nlon, isym, mtrunc, x, y, idxy, wshp, lwshp, iwshp, liwshp, work, lwork, ierror)
    !
    !     shpg projects the array x onto the set of functions represented
    !     by a discrete set of spherical harmonics.
    !
    !     input parameters
    !
    !     nlat   the number of colatitudes on the full sphere including the
    !            poles. for example, nlat = 37 for a five degree grid.
    !            nlat determines the grid increment in colatitude as
    !            pi/(nlat-1).  if nlat is odd the equator is located at
    !            grid point i=(nlat + 1)/2. if nlat is even the equator is
    !            located half way between points i=nlat/2 and i=nlat/2+1.
    !            nlat must be at least 3.
    !
    !     nlon   the number of distinct londitude points.  nlon determines
    !            the grid increment in longitude as 2*pi/nlon. for example
    !            nlon = 72 for a five degree grid. nlon must be greater
    !            than or equal to 4. the efficiency of the computation is
    !            improved when nlon is a product of small prime numbers.
    !            nlon must be at least 4.
    !
    !     isym   currently not used.
    !
    !     mtrunc the highest longitudinal wave number retained in the
    !            projection. It must be less than or equal to
    !            the minimum of nlat-1 and nlon/2. The first wave
    !            number is zero. For example, if wave numbers 0 and
    !            1 are desired then mtrunc = 1.

    !            zero.
    !
    !     x      a two dimensional array that contains the the nlat
    !            by nlon array x(i, j) defined at the colatitude point
    !            theta(i) = (i-1)*pi/(nlat-1) and longitude point phi(j) =
    !            (j-1)*2*pi/nlon.
    !
    !     idxy   the first dimension of the arrays x and y as they
    !            appear in the program that calls shpg. It must be
    !            at least nlat.
    !
    !     wshp   a single precision array that must be saved for
    !            repeated use by subroutine shpg.
    !
    !     lwshp  the dimension of the array wshp as it appears in the
    !            program that calls shpgi. It must be at least
    !            2*(nlat + 1)**2+nlon+log2(nlon)
    !
    !     iwshp  an integer array that must be saved for repeated
    !            use by subroutine shpg.
    !
    !
    !     liwshp the dimension of the array iwshp as it appears in the
    !            program that calls shpgi. It must be at least
    !            4*(nlat + 1).
    !
    !     work   a single precision work array that does
    !            not have to be saved.
    !
    !     lwork  the dimension of the array work as it appears in the
    !            program that calls shpg. It must be at least
    !            max(nlat*nlon, 4*(nlat + 1)).
    !
    !     **************************************************************
    !
    !     output parameters
    !
    !     y      an nlat by nlon single precision array that contains
    !            the projection of x onto the set of functions that
    !            can be represented by the discrete set of spherical
    !            harmonics. The arrays x(i, j) and y(i, j) are located
    !            at colatitude point theta(i) = (i-1)*pi/(nlat-1) and
    !            longitude point phi(j) = (j-1)*2*pi/nlon.
    !
    !     ierror = 0  no errors
    !            = 1  error in the specification of nlat
    !            = 2  error in the specification of nlon
    !            = 3  error in the specification of isym
    !            = 4  error in the specification of mtrunc
    !            = 5  error in the specification of lwshp
    !            = 6  error in the specification of liwshp
    !            = 7  error in the specification of lwork
    !
    module subroutine shpg(nlat, nlon, isym, mtrunc, x, y, idxy, &
        wshp, lwshp, iwshp, liwshp, work, lwork, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: isym
        integer(ip), intent(in)  :: mtrunc
        real(wp),    intent(in)  :: x(idxy, nlon)
        real(wp),    intent(out) :: y(idxy, nlon)
        integer(ip), intent(in)  :: idxy
        real(wp),    intent(in)  :: wshp(lwshp)
        integer(ip), intent(in)  :: lwshp
        integer(ip), intent(in)  :: iwshp(liwshp)
        integer(ip), intent(in)  :: liwshp
        real(wp),    intent(out) :: work(lwork)
        integer(ip), intent(in)  :: lwork
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip) :: iw1
        integer(ip) :: iw2
        integer(ip) :: iw3
        integer(ip) :: iw4
        integer(ip) :: jw1
        integer(ip) :: jw2
        integer(ip) :: jw3
        integer(ip) :: jw4
        integer(ip) :: log2n
        integer(ip) :: lw1
        integer(ip) :: mmax
        integer(ip) :: mwrk
        integer(ip) :: nloc1, nte
        integer(ip) :: nloc2
        type(SpherepackUtility) :: util

        ! Check calling arguments
        ierror = 1
        if (nlat < 1) return
        ierror = 2
        if (nlon < 1) return
        ierror = 3
        if (isym < 0 .or. isym > 2) return
        ierror = 4
        mmax = min(nlat-1, nlon/2)
        if (mtrunc<0 .or. mtrunc>mmax) return
        ierror = 5
        log2n = int(log(real(nlon, kind=wp))/log(TWO), kind=ip)
        lw1 = 2*(nlat + 1)**2
        if (lwshp<lw1+nlon+log2n) return
        ierror = 6
        if (liwshp<4*(nlat + 1)) return
        ierror = 7
        mwrk = max(nlat*nlon, 4*(nlat + 1))
        if (lwork <mwrk) return
        ierror = 0

        y(1:nlat, :) = x(1:nlat, :)

        call util%hfft%forward(nlat, nlon, y, idxy, wshp(lw1+1))

        ! Set workspace index pointers
        nte = (nlat + 1)/2
        nloc1 = 2*(nte**2)
        nloc2 = nlat+1
        lw1 = 2*(nlat + 1)**2
        iw1 = 1
        iw2 = iw1+nloc1
        iw3 = iw2+nloc1
        iw4 = iw3+nloc1
        jw1 = 1
        jw2 = jw1+nloc2
        jw3 = jw2+nloc2
        jw4 = jw3+nloc2

        call shpg_lower_utility_routine(nlat, nlon, isym, mtrunc, y, y, idxy, ierror, &
            nte, wshp(iw1), wshp(iw2), wshp(iw3), wshp(iw4), iwshp(jw1), &
            iwshp(jw2), iwshp(jw3), iwshp(jw4), work(jw1), &
            work(jw2), work(jw3), work(jw4))

        call util%hfft%backward(nlat, nlon, y, idxy, wshp(lw1+1))

        y(1: nlat, :) = y(1:nlat, :)/nlon

    end subroutine shpg

    ! Purpose:
    !
    !     shpgi initializes the arrays wshp and iwshp for subsequent
    !     use in subroutine shpg, which performs the harmonic projection
    !     which is equivalent to a harmonic analysis followed by
    !     harmonic synthesis but faster and with less memory.
    !     (see description of subroutine shpg below).
    !
    !     subroutine shpgi(nlat, nlon, isym, mtrunc, wshp, lwshp, iwshp, liwshp, work, lwork, ierror)
    !
    !     shpgi initializes arrays wshp and iwshp for repeated use
    !     by subroutine shpg ....
    !
    !     input parameters
    !
    !     nlat   the number of colatitudes on the full sphere including the
    !            poles. for example, nlat = 37 for a five degree grid.
    !            nlat determines the grid increment in colatitude as
    !            pi/(nlat-1).  if nlat is odd the equator is located at
    !            grid point i=(nlat + 1)/2. if nlat is even the equator is
    !            located half way between points i=nlat/2 and i=nlat/2+1.
    !            nlat must be at least 3.
    !
    !     nlon   the number of distinct londitude points.  nlon determines
    !            the grid increment in longitude as 2*pi/nlon. for example
    !            nlon = 72 for a five degree grid. nlon must be greater
    !            than or equal to 4. the efficiency of the computation is
    !            improved when nlon is a product of small prime numbers.
    !            nlon must be at least 4.
    !
    !     isym   currently not used, no equatorial symmetries assumed, 
    !            only whole sphere computations.
    !
    !     mtrunc the highest longitudinal wave number retained in the
    !            projection. It must be less than or equal to
    !            the minimum of nlat-1 and nlon/2. The first wave
    !            number is zero. For example, if wave numbers 0 and
    !            1 are desired then mtrunc = 1.
    !
    !     lwshp  the dimension of the array wshp as it appears in the
    !            program that calls shpgi. It must be at least
    !            2*(nlat + 1)**2+nlon+log2(nlon)
    !
    !     liwshp the dimension of the array iwshp as it appears in the
    !            program that calls shpgi. It must be at least
    !            4*(nlat + 1).
    !
    !     lwork  the dimension of the array work as it appears in the
    !            program that calls shpgi. It must be at least
    !            1.25*(nlat + 1)**2+7*nlat+8.
    !
    !     **************************************************************
    !
    !     output parameters
    !
    !     wshp   a single precision array that must be saved for
    !            repeated use by subroutine shpg.
    !
    !     iwshp  an integer array that must be saved for repeated
    !            use by subroutine shpg.
    !
    !     work   a real work array that does
    !            not have to be saved.
    !
    !     ierror = 0  no errors
    !            = 1  error in the specification of nlat
    !            = 2  error in the specification of nlon
    !            = 3  error in the specification of isym
    !            = 4  error in the specification of mtrunc
    !            = 5  error in the specification of lwshp
    !            = 6  error in the specification of liwshp
    !            = 7  error in the specification of lwork
    !
    module subroutine shpgi(nlat, nlon, isym, mtrunc, wshp, lwshp, iwshp, &
        liwshp, work, lwork, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: isym
        integer(ip), intent(in)  :: mtrunc
        real(wp),    intent(out) :: wshp(lwshp)
        integer(ip), intent(in)  :: lwshp
        integer(ip), intent(in)  :: iwshp(liwshp)
        integer(ip), intent(in)  :: liwshp
        real(wp),    intent(out) :: work(lwork)
        integer(ip), intent(in)  :: lwork
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip) :: iw1
        integer(ip) :: iw2
        integer(ip) :: iw3
        integer(ip) :: iw4
        integer(ip) :: jw1
        integer(ip) :: jw2
        integer(ip) :: jw3
        integer(ip) :: jw4
        integer(ip) :: ktot
        integer(ip) :: kw1
        integer(ip) :: kw10
        integer(ip) :: kw11
        integer(ip) :: kw2
        integer(ip) :: kw3
        integer(ip) :: kw4
        integer(ip) :: kw5
        integer(ip) :: kw6
        integer(ip) :: kw7
        integer(ip) :: kw8
        integer(ip) :: kw9
        integer(ip) :: log2n
        integer(ip) :: lw1
        integer(ip) :: mlwk
        integer(ip) :: mmax
        integer(ip) :: nloc1
        integer(ip) :: nloc2
        integer(ip) :: nte
        type(SpherepackUtility) :: util

        ! Check calling arguments
        ierror = 1
        if (nlat<1) return
        ierror = 2
        if (nlon<1) return
        ierror = 3
        if (isym < 0 .or. isym > 2) return
        ierror = 4
        mmax = min(nlat-1, nlon/2)
        if (mtrunc<0 .or. mtrunc>mmax) return
        ierror = 5
        lw1 = 2*((nlat + 1)**2)
        log2n = int(log(real(nlon, kind=wp))/log(TWO), kind=ip)
        if (lwshp<lw1+nlon+log2n) return
        ierror = 6
        if (liwshp<4*(nlat + 1)) return
        ierror = 7
        mlwk = (5*((nlat + 1)**2+7*nlat+8))/4
        if (lwork <mlwk) return
        ierror = 0

        call util%hfft%initialize(nlon, wshp(lw1+1))

        nte = (nlat + 1)/2
        nloc1 = 2*nte*nte
        nloc2 = nlat+1
        lw1 = 2*((nlat + 1)**2)
        iw1 = 1
        iw2 = iw1+nloc1
        iw3 = iw2+nloc1
        iw4 = iw3+nloc1
        jw1 = 1
        jw2 = jw1+nloc2
        jw3 = jw2+nloc2
        jw4 = jw3+nloc2
        kw1 = 1
        kw2 = kw1+nte
        kw3 = kw2+nte
        kw4 = kw3+2*nte
        kw5 = kw4+2*nte
        kw6 = kw5+nte
        kw7 = kw6+nte
        kw8 = kw7+4*nte
        kw9 = kw8+2*nte
        kw10 = kw9+nloc1
        kw11 = kw10+nloc1
        ktot = kw11+nte*nte

        call shpgi_lower_utility_routine(nlat, nlon, isym, mtrunc, nte, ierror, wshp(iw1), wshp(iw2), &
            wshp(iw3), wshp(iw4), iwshp(jw1), iwshp(jw2), iwshp(jw3), &
            iwshp(jw4), work(kw1), work(kw2), work(kw3), work(kw4), work(kw5), &
            work(kw6), work(kw7), work(kw8), work(kw9), work(kw10), work(kw11))

    end subroutine shpgi

    subroutine shpgi_lower_utility_routine(nlat, nlon, isym, mtrunc, idp, ierror, &
        pe, po, ze, zo, ipse, jzse, ipso, jzso, &
        cp, wx, thet, gwts, xx, z, a, b, ped, pod, u)

        real(wp) :: dfn
        real(wp) :: dmax
        integer(ip) :: i
        integer(ip) :: idp
        integer(ip) :: ierr
        integer(ip) :: ierror
        integer(ip) :: iip
        integer(ip) :: ipse
        integer(ip) :: ipso
        integer(ip) :: isym
        integer(ip) :: it
        integer(ip) :: j
        integer(ip) :: js
        integer(ip) :: jzse
        integer(ip) :: jzso
        integer(ip) :: k
        integer(ip) :: lock
        
        integer(ip) :: m
        integer(ip) :: modn
        integer(ip) :: mp1
        integer(ip) :: ms2
        integer(ip) :: mtrunc
        integer(ip) :: mxtr
        integer(ip) :: n
        integer(ip) :: nec
        integer(ip) :: nem
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: nmx
        integer(ip) :: noc
        integer(ip) :: nom
        integer(ip) :: ns2
        integer(ip) :: nshe
        integer(ip) :: nsho
        integer(ip) :: nte
        integer(ip) :: nto
        real(wp) :: pe
        real(wp) :: po
        real(wp) :: sum1
        real(wp) :: toe
        real(wp) :: tusl
        real(wp) :: ze
        real(wp) :: zo
        real(wp) :: zort

        real(wp) :: summation, a1, b1, c1
        real(wp) :: cp(idp), wx(idp), &
            thet(nlat), gwts(nlat), xx(idp), z(idp), a(4*idp), &
            b(2*idp), ped(idp, idp, 2), pod(idp, idp, 2), u(idp, idp)

        dimension pe(idp, idp, 2), po(idp, idp, 2), ze(idp, idp, 2), &
            zo(idp, idp, 2), &
            ipse(idp, 2), jzse(idp, 2), ipso(idp, 2), jzso(idp, 2), &
            nshe(2), nsho(2)
        dimension zort(64, 64, 2)

        type(SpherepackUtility) :: util

        ns2 = nlat/2
        modn = nlat-2*ns2
        nte = (nlat + 1)/2
        nto = nlat-nte
        tusl = ZERO
        toe = ZERO

        ! Compute gauss grid distribution
        call compute_gaussian_latitudes_and_weights(nlat, thet, gwts, ierr)

        gwts(1:nto) = TWO * gwts(1:nto)

        ! Compute n**2 basis (even functions)
        do n=1, 2*nlat-2
            dfn = n
            a(n) = sqrt(dfn * (dfn + ONE))
        end do

        do n=1, nlat-1
            dfn = n
            b(n) = sqrt((TWO * dfn + THREE)/(TWO * dfn - ONE))
        end do

        mxtr = min(nlat-1, nlon/2, mtrunc)
        iip = 2
        generate_even_functions: do mp1=1, mxtr+1
            m = mp1-1
            iip = 3-iip
            ms2 = mp1/2
            nem = (nlat-m+1)/2
            nec = nte-nem

            ! Compute associated legendre functions
            if (m <= 1) then
                do j=1, nem
                    n = 2*j+m-2
                    call util%compute_fourier_coefficients(m, n, cp)
                    do i=1, nte
                        call util%compute_legendre_polys_from_fourier_coeff(m, n, thet(i), cp, ped(i, j+nec, iip))
                    end do
                end do
            else
                do j=1, nem
                    n = 2*j+m-2
                    if (m>1 .and. n>mxtr) then
                        do i=1, nte
                            u(i, j+nec) = ped(i, j+nec, iip)
                        end do
                    else
                        a1 = b(n-1)*a(n+m-3)/a(n+m-1)
                        b1 = a(n-m+1)/a(n+m-1)
                        if (n-m<=1) then
                            do i=1, nte
                                u(i, j+nec) = a1*ped(i, j+nec-1, iip) &
                                    - b1*ped(i, j+nec, iip)
                            end do
                        else
                            c1 = b(n-1)*a(n-m-1)/a(n+m-1)
                            do i=1, nte
                                u(i, j+nec) = a1*ped(i, j+nec-1, iip) &
                                    - b1*ped(i, j+nec, iip) + c1*u(i, j+nec-1)
                            end do
                        end if
                    end if
                end do
                do j=1, nem
                    do i=1, nte
                        ped(i, j+nec, iip) = u(i, j+nec)
                    end do
                end do
            end if
            if (nec<=0) cycle generate_even_functions
            !
            !     generate orthogonal vector with
            !     random numbers
            call random_seed()
            call random_number(xx(1:nte))

            it = 0
            generate_random_orth_vec_even: do
                it = it+1
                if (it > 2) exit generate_random_orth_vec_even

                do i=1, nte
                    z(i) = ZERO
                    wx(i) = gwts(i)*xx(i)
                end do

                do j=1, nte
                    if (j == nec) cycle
                    call accumulate_inner_products(nte, wx, ped(1, j, iip), z)
                end do

                do i=1, nte
                    xx(i) = xx(i)-z(i)
                end do

                call calculate_normal(nte, xx, idp, gwts)

            end do generate_random_orth_vec_even

            do i=1, nte
                ped(i, nec, iip) = xx(i)
            end do
        end do generate_even_functions
        !
        !     reorder if mtrunc is less than nlat-1
        !         case of even functions
        !
        nmx = nlat-mxtr
        if (modn==1) then
            nshe(1) = nmx/2
            nshe(2) = (nmx-1)/2
        else
            nshe(1) = (nmx-1)/2
            nshe(2) = nmx/2
        end if
        !
        do mp1=1, 2
            do j=1, nte
                js = j+nshe(mp1)
                if (js>nte) js = js-nte
                do i=1, nte
                    u(i, js) = ped(i, j, mp1)
                end do
            end do
            do j=1, nte
                do i=1, nte
                    ped(i, j, mp1) = u(i, j)
                end do
            end do
        end do
        !
        call truncate(0, nte, idp, ped(1, 1, 1), nte, ipse(1, 1))
        call truncate(0, nte, idp, ped(1, 1, 2), nte, ipse(1, 2))
        !
        ! Compute the analysis matrices
        !
        do iip=1, 2
            do i=1, nte
                lock = 0
                do j=1, nte
                    summation = ped(j, i, iip)*gwts(j)
                    ze(j, i, iip) =  summation
                    pe(i, j, iip) = ped(i, j, iip)
                    if (abs(summation)>MACHINE_EPSILON .and. lock==0) then
                        lock = 1
                        jzse(i, iip) = j
                    end if
                end do
            end do
        end do
        !
        !     check orthogonality of pe(i, j, mp1)  mp1=1, 2
        !
        do iip=1, 2
            dmax = ZERO
            do i=1, nte
                do j=1, nte
                    sum1 = ZERO
                    do k=1, nte
                        sum1 = sum1+ze(k, i, iip)*pe(k, j, iip)
                    end do
                    zo(i, j, iip) = sum1
                    if (i/=j) then
                        dmax = max(dmax, abs(sum1))
                    else
                        dmax = max(dmax, abs(sum1-ONE))
                    end if
                end do
            end do
        end do
        !
        ! Compute n**2 basis (odd functions)
        !
        iip = 2
        generate_odd_functions: do mp1=1, mxtr+1
            iip = 3-iip
            m = mp1-1
            ms2 = mp1/2
            nem = (nlat-m+1)/2
            nom = nlat-m-nem
            noc = nto-nom
            !
            ! Compute associated legendre functions
            !
            if (m<=1) then
                do j=1, nom
                    n = 2*j+m-1
                    call util%compute_fourier_coefficients(m, n, cp)
                    do i=1, nte
                        call util%compute_legendre_polys_from_fourier_coeff(m, n, thet(i), cp, pod(i, j+noc, iip))
                    end do
                    if (modn>0) pod(nte, j+noc, iip) = ZERO
                end do
            else
                do j=1, nom
                    n = 2*j+m-1
                    if (m>1 .and. n>mxtr) then
                        do i=1, nte
                            u(i, j+noc) = pod(i, j+noc, iip)
                        end do
                    else
                        a1 = b(n-1)*a(n+m-3)/a(n+m-1)
                        b1 = a(n-m+1)/a(n+m-1)
                        if (n-m<=1) then
                            do i=1, nte
                                u(i, j+noc) = a1*pod(i, j+noc-1, iip) &
                                    - b1*pod(i, j+noc, iip)
                            end do
                        else
                            c1 = b(n-1)*a(n-m-1)/a(n+m-1)
                            do i=1, nte
                                u(i, j+noc) = a1*pod(i, j+noc-1, iip) &
                                    - b1*pod(i, j+noc, iip) + c1*u(i, j+noc-1)
                            end do
                        end if
                    end if
                    if (modn==1) u(nte, j+noc) = ZERO
                end do
                do j=1, nom
                    do i=1, nte
                        pod(i, j+noc, iip) = u(i, j+noc)
                    end do
                end do
            end if

            if (noc <= 0) cycle generate_odd_functions

            call random_number(xx(1:nte))

            if (modn==1) xx(nte) = ZERO

            it = 0
            generate_random_orth_vec_odd: do
                it = it+1
                if (it > 2) exit generate_random_orth_vec_odd

                z(1: nte) = ZERO
                wx(1: nte) = gwts(1: nte)*xx(1: nte)

                do j=1, nto
                    if (j==noc) cycle
                    call accumulate_inner_products(nte, wx, pod(1, j, iip), z(1))
                end do

                xx(1: nte) = xx(1: nte)-z(1: nte)

                call calculate_normal(nte, xx, idp, gwts)
            end do generate_random_orth_vec_odd

            pod(1: nte, noc, iip) = xx(1: nte)

            if (modn==1) pod(nte, noc, iip) = ZERO
        end do generate_odd_functions

        nmx = nlat-mxtr

        if (modn==1) then
            nsho(1) = (nmx-1)/2
            nsho(2) = nmx/2
        else
            nsho(1) = nmx/2
            nsho(2) = (nmx-1)/2
        end if
        !
        do mp1=1, 2
            do j=1, nto
                js = j+nsho(mp1)
                if (js>nto) js = js-nto
                do i=1, nte
                    u(i, js) = pod(i, j, mp1)
                end do
            end do
            do j=1, nto
                do i=1, nte
                    pod(i, j, mp1) = u(i, j)
                end do
            end do
        end do
        !
        call truncate(0, nte, idp, pod(1, 1, 1), nto, ipso(1, 1))
        call truncate(0, nte, idp, pod(1, 1, 2), nto, ipso(1, 2))
        !
        ! Compute the analysis matrices (odd functions)
        !
        do iip=1, 2
            do i=1, nto
                lock = 0
                do j=1, nto
                    summation = pod(j, i, iip)*gwts(j)
                    zo(j, i, iip) = summation
                    po(i, j, iip) = pod(i, j, iip)
                    if (abs(summation)>MACHINE_EPSILON .and. lock==0) then
                        lock = 1
                        jzso(i, iip) = j
                    end if
                end do
            end do
        end do
        !
        !     check orthogonality of po(i, j, mp1)  mp1=1, 2
        !
        do iip=1, 2
            dmax = ZERO
            do i=1, nto
                do j=1, nto
                    sum1 = ZERO
                    do k=1, nto
                        sum1 = sum1+zo(k, i, iip)*po(k, j, iip)
                    end do
                    zort(i, j, iip) = sum1
                    if (i/=j) then
                        dmax = max(dmax, abs(sum1))
                    else
                        dmax = max(dmax, abs(sum1-ONE))
                    end if
                end do
            end do
        end do

    end subroutine shpgi_lower_utility_routine

    subroutine shpg_lower_utility_routine(nlat, nlon, isym, mtrunc, sx, sy, idxy, ierror, &
        idp, pe, po, ze, zo, ipse, jzse, ipso, jzso, xe, xo, ye, yo)

        integer(ip) :: i
        integer(ip) :: idp
        integer(ip) :: idxy
        integer(ip) :: ierror
        integer(ip) :: iip
        integer(ip) :: ipse
        integer(ip) :: ipso
        integer(ip) :: isym
        integer(ip) :: j
        integer(ip) :: js
        integer(ip) :: jzse
        integer(ip) :: jzso
        integer(ip) :: lag
        integer(ip) :: m
        integer(ip) :: modn
        integer(ip) :: mp1
        integer(ip) :: mpm
        integer(ip) :: ms2
        integer(ip) :: mtrunc
        integer(ip) :: mxtr
        integer(ip) :: nec
        integer(ip) :: nem
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: nmx
        integer(ip) :: noc
        integer(ip) :: nom
        integer(ip) :: ns2
        integer(ip) :: nshe
        integer(ip) :: nsho
        integer(ip) :: nte
        integer(ip) :: nto
        real(wp) :: pe
        real(wp) :: po
        real(wp) :: sx
        real(wp) :: sy
        real(wp) :: xe
        real(wp) :: xo
        real(wp) :: ye
        real(wp) :: yo
        real(wp) :: ze
        real(wp) :: zo
        !
        dimension sx(idxy, nlon), sy(idxy, nlon), nshe(2), nsho(2), &
            pe(idp, idp, 2), po(idp, idp, 2), ze(idp, idp, 2), zo(idp, idp, 2), &
            ipse(idp, 2), jzse(idp, 2), ipso(idp, 2), jzso(idp, 2), &
            xe(idp, 2), xo(idp, 2), ye(idp, 2), yo(idp, 2)
        !
        ns2 = nlat/2
        modn = nlat-ns2-ns2
        nte = (nlat + 1)/2
        nto = nlat-nte
        !
        mxtr = min(nlat-1, nlon/2, mtrunc)
        nmx = nlat-mxtr
        if (modn==1) then
            nshe(1) = nmx/2
            nshe(2) = (nmx-1)/2
            nsho(1) = (nmx-1)/2
            nsho(2) = nmx/2
        else
            nshe(1) = (nmx-1)/2
            nshe(2) = nmx/2
            nsho(1) = nmx/2
            nsho(2) = (nmx-1)/2
        end if
        !
        iip = 2
        outer_loop: do mp1=1, mxtr+1
            iip = 3-iip
            if (mxtr==nlat-1 .and. mp1==1) then
                do i=1, nlat
                    sy(i, mp1) = sx(i, mp1)
                end do
                cycle outer_loop
            end if
            m = mp1-1
            mpm = max(1, m+m)
            ms2 = mp1/2
            nem = (nlat-m+1)/2-nshe(iip)
            nom = (nlat-m)/2-nsho(iip)
            nec = nte-nem
            noc = nto-nom

            do i=1, nte
                xe(i, 1) = HALF * (sx(i, mpm)+sx(nlat+1-i, mpm))
                xo(i, 1) = HALF * (sx(i, mpm)-sx(nlat+1-i, mpm))
            end do

            if (mpm<nlon) then
                do i=1, nte
                    xe(i, 2) = HALF * (sx(i, mpm+1)+sx(nlat+1-i, mpm+1))
                    xo(i, 2) = HALF * (sx(i, mpm+1)-sx(nlat+1-i, mpm+1))
                end do
            end if

            lag = 0
            if (m==0.or.mpm==nlon) lag = 1
            if (3*nec<2*nem.or.nem==0) then
                call matrix_multiply(lag, nte, nec, idp, pe(1, 1, iip), nte, idp, &
                    ze(1, 1, iip), xe, ye, ipse(1, iip), jzse(1, iip))
                do i=1, nte
                    ye(i, 1) = xe(i, 1)-ye(i, 1)
                end do
                if (mpm<nlon .and. m/=0) then
                    do i=1, nte
                        ye(i, 2) = xe(i, 2)-ye(i, 2)
                    end do
                end if
            else
                call matrix_multiply(lag, nte, nem, idp, pe(1, nec+1, iip), nte, idp, &
                    ze(1, nec+1, iip), xe, ye, ipse(nec+1, iip), jzse(nec+1, iip))
            end if
            if (3*noc<2*nom.or.nom==0) then
                call matrix_multiply(lag, nto, noc, idp, po(1, 1, iip), nto, idp, &
                    zo(1, 1, iip), xo, yo, ipso(1, iip), jzso(1, iip))
                do i=1, nto
                    yo(i, 1) = xo(i, 1)-yo(i, 1)
                end do
                if (mpm<nlon .and. m/=0) then
                    do i=1, nto
                        yo(i, 2) = xo(i, 2)-yo(i, 2)
                    end do
                end if
            else
                call matrix_multiply(lag, nto, nom, idp, po(1, noc+1, iip), nto, idp, &
                    zo(1, noc+1, iip), xo, yo, ipso(noc+1, iip), jzso(noc+1, iip))
            end if
            do i=1, nto
                sy(i, mpm) = ye(i, 1)+yo(i, 1)
                sy(nlat+1-i, mpm) = ye(i, 1)-yo(i, 1)
            end do
            if (nte>nto) sy(nte, mpm) = ye(nte, 1)
            if (mpm<nlon .and. m/=0) then
                do i=1, nto
                    sy(i, mpm+1) = ye(i, 2)+yo(i, 2)
                    sy(nlat+1-i, mpm+1) = ye(i, 2)-yo(i, 2)
                end do
                if (nte>nto) sy(nte, mpm+1) = ye(nte, 2)
            end if
        end do outer_loop

        js = mxtr+mxtr+2
        do j=js, nlon
            do i=1, nlat
                sy(i, j) = ZERO
            end do
        end do

    end subroutine shpg_lower_utility_routine

    subroutine matrix_multiply(lag, lr, lc, ld, a, mc, md, b, x, y, is, js)

        ! Dummy arguments
        integer(ip), intent(in)  :: lag
        integer(ip), intent(in)  :: lc
        integer(ip), intent(in)  :: ld
        integer(ip), intent(in)  :: lr
        real(wp),    intent(in)  :: a(ld, *)
        integer(ip), intent(in)  :: mc
        integer(ip), intent(in)  :: md
        real(wp),    intent(in)  :: b(md, *)
        real(wp),    intent(in)  :: x(ld, 2)
        real(wp),    intent(inout) :: y(ld, 2)
        integer(ip), intent(in)  :: is(*)
        integer(ip), intent(in)  :: js(*)

        ! Local variables
        integer(ip) :: i, j, k, kmx
        real(wp)    :: sum1, sum2

        kmx = min(lr+1, ld)

        select case (lag)
            case(1)
                y(1: kmx, 1) = ZERO
                if (lc > 0) then
                    do i=1, lc
                        sum1 = ZERO
                        do j=js(i), mc
                            sum1 = sum1 + b(j, i)*x(j, 1)
                        end do
                        do k=is(i), lr
                            y(k, 1) = y(k, 1)+sum1*a(k, i)
                        end do
                    end do
                end if
            case default
                y(1: kmx, :) = ZERO
                if (lc > 0) then
                    do i=1, lc
                        sum1 = ZERO
                        sum2 = ZERO
                        do j=js(i), mc
                            sum1 = sum1 + b(j, i)*x(j, 1)
                            sum2 = sum2 + b(j, i)*x(j, 2)
                        end do
                        do k=is(i), lr
                            y(k, 1) = y(k, 1)+sum1*a(k, i)
                            y(k, 2) = y(k, 2)+sum2*a(k, i)
                        end do
                    end do
                end if
        end select

    end subroutine matrix_multiply

    subroutine calculate_normal(n, x, id, q)

        ! Dummy arguments
        integer(ip), intent(in)    :: n
        real(wp),    intent(inout) :: x(n)
        integer(ip), intent(in)    :: id
        real(wp),    intent(in)    :: q(n)

        ! Local variables
        integer(ip) :: i
        real(wp)    :: sqs

        ! Normalize x
        sqs = ZERO
        do i=1, n
            sqs = sqs+q(i)*(x(i)**2)
        end do

        x = x/sqrt(sqs)

    end subroutine calculate_normal

end submodule scalar_projection_gaussian_grid
