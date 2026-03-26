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
!     This file contains code and documentation for subroutines
!     shpei and shpe.
!
submodule(scalar_projection_routines) scalar_projection_regular_grid

contains

    ! Purpose:
    !
    !     the n**2 projection with complement, odd/even
    !     factorization and zero truncation on an
    !     equally spaced grid as defined in the JCP paper
    !     "Generalized discrete spherical harmonic transforms"
    !     by Paul N. Swarztrauber and William F. Spotz
    !     It is equivalent to a harmonic analysis followed
    !     by a synthesis except faster and requires less memory.
    !
    !     subroutine shpe(nlat, nlon, isym, mtrunc, x, y, idxy, &
    !     wshp, lwshp, iwshp, liwshp, work, lwork, ierror)
    !
    !     shpe projects the array x onto the set of functions represented
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
    !            nlon must beat least 4.
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
    !            appear in the program that calls shpe. It must be
    !            at least nlat.
    !
    !     wshp   a single precision array that must be saved for
    !            repeated use by subroutine shpe.
    !
    !     lwshp  the dimension of the array wshp as it appears in the
    !            program that calls shpei. It must be at least
    !            2*(nlat + 1)**2+nlon+log2(nlon)
    !
    !     iwshp  an integer array that must be saved for repeated
    !            use by subroutine shpe.
    !
    !
    !     liwshp the dimension of the array iwshp as it appears in the
    !            program that calls shpei. It must be at least
    !            4*(nlat + 1).
    !
    !     work   a single precision work array that does
    !            not have to be saved.
    !
    !     lwork  the dimension of the array work as it appears in the
    !            program that calls shpe. It must be at least
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
    module subroutine shpe(nlat, nlon, isym, mtrunc, x, y, idxy, &
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
        integer(ip) :: nloc1
        integer(ip) :: nloc2
        integer(ip) :: nte
        type(SpherepackUtility) :: util

        ! Check calling arguments
        ierror = 1
        if (nlat < 3) return
        ierror = 2
        if (nlon < 4) return
        ierror = 3
        if (isym < 0 .or. isym > 2) return
        ierror = 4
        mmax = min(nlat-1, nlon/2)
        if (mtrunc < 0 .or. mtrunc > mmax) return
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
        nloc1 = 2*nte*nte
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

        call shpe_lower_utility_routine(nlat, nlon, isym, mtrunc, y, y, idxy, ierror, &
            nte, wshp(iw1), wshp(iw2), wshp(iw3), wshp(iw4), iwshp(jw1), &
            iwshp(jw2), iwshp(jw3), iwshp(jw4), work(jw1), &
            work(jw2), work(jw3), work(jw4))

        call util%hfft%backward(nlat, nlon, y, idxy, wshp(lw1+1))

        y(1:nlat, :) = y(1:nlat, :)/nlon

    end subroutine shpe

    ! Purpose:
    !
    !     Initializes arrays wshp and iwshp for
    !     subsequent repeated use by subroutine shpe, which
    !     performs the harmonic projection equivalent to a
    !     harmonic analysis followed by harmonic synthesis
    !     but faster and with less memory.
    !
    !     subroutine shpei(nlat, nlon, isym, mtrunc, wshp, lwshp, iwshp, 
    !     liwshp, work, lwork, ierror)
    !
    !     shpei initializes arrays wshp and iwshp for repeated use
    !     by subroutine shpe ....
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
    !            nlon must beat least 4.
    !
    !     isym   currently not used.
    !
    !     mtrunc the highest longitudinal wave number retained in the
    !            projection. It must be less than or equal to
    !            the minimum of nlat-1 and nlon/2. The first wave
    !            number is zero. For example, if wave numbers 0 and
    !            1 are desired then mtrunc = 1.
    !
    !     lwshp  the dimension of the array wshp as it appears in the
    !            program that calls shpei. It must be at least
    !            2*(nlat + 1)**2+nlon+log2(nlon)
    !
    !     liwshp the dimension of the array iwshp as it appears in the
    !            program that calls shpei. It must be at least
    !            4*(nlat + 1).
    !
    !     lwork  the dimension of the array work as it appears in the
    !            program that calls shpei. It must be at least
    !            1.25*(nlat + 1)**2+7*nlat+8.
    !
    !     **************************************************************
    !
    !     output parameters
    !
    !     wshp   a single precision array that must be saved for
    !            repeated use by subroutine shpe.
    !
    !     iwshp  an integer array that must be saved for repeated
    !            use by subroutine shpe.
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
    module subroutine shpei(nlat, nlon, isym, mtrunc, wshp, lwshp, iwshp, &
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
        integer(ip) :: kw1
        integer(ip) :: kw10
        integer(ip) :: kw11
        integer(ip) :: kw12
        integer(ip) :: kw13
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

        ! Set contants
        mmax = min(nlat-1, nlon/2)
        lw1 = 2*(nlat + 1)**2
        log2n = int(log(real(nlon, kind=wp))/log(TWO), kind=ip)
        mlwk = 5*((nlat + 1)**2 + 7*nlat + 8)/4

        ! Check calling arguments
        if (nlat < 3) then
            ierror = 1
        else if (nlon < 4) then
            ierror = 2
        else if (isym < 0 .or. isym > 2) then
            ierror = 3
        else if (mtrunc < 0 .or. mtrunc > mmax) then
            ierror = 4
        else if (lwshp < lw1+nlon+log2n) then
            ierror = 5
        else if (liwshp < 4*(nlat + 1)) then
            ierror = 6
        else if (lwork < mlwk) then
            ierror = 7
        else
            ierror = 0
        end if

        ! Check error flag
        if (ierror /= 0) return

        call util%hfft%initialize(nlon, wshp(lw1+1))

        ! Set workspace index pointers
        nte = (nlat + 1)/2
        nloc1 = 2*nte*nte
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
        kw1 = 1
        kw2 = kw1+nte
        kw3 = kw2+nte
        kw4 = kw3+nte
        kw5 = kw4+nte+1
        kw6 = kw5+nte
        kw7 = kw6+nte
        kw8 = kw7+nte
        kw9 = kw8+nte
        kw10 = kw9+nloc2+nloc2
        kw11 = kw10+nloc2
        kw12 = kw11+nloc1
        kw13 = kw12+nloc1

        call shpei_lower_utility_routine(nlat, nlon, isym, mtrunc, nte, ierror, wshp(iw1), wshp(iw2), &
            wshp(iw3), wshp(iw4), iwshp(jw1), iwshp(jw2), iwshp(jw3), &
            iwshp(jw4), work(kw1), work(kw2), work(kw3), work(kw4), work(kw5), &
            work(kw6), work(kw7), work(kw8), work(kw9), work(kw10), work(kw11), &
            work(kw12), work(kw11), work(kw12), work(kw13))

    end subroutine shpei

    subroutine shpei_lower_utility_routine(nlat, nlon, isym, mtrunc, idp, ierror, &
        pe, po, ze, zo, ipse, jzse, ipso, jzso, &
        cp, work, wx, s, e, thet, xx, z, a, b, we, ped, wo, pod, u)

        real(wp) :: dfn
        integer(ip) :: i
        integer(ip) :: idp
        integer(ip) :: ierror
        integer(ip) :: info
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
        integer(ip) :: mrank
        integer(ip) :: ms2
        integer(ip) :: mtrunc
        integer(ip) :: mxtr
        integer(ip) :: n
        integer(ip) :: nem
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: nom
        integer(ip) :: nrank
        integer(ip) :: ns2
        integer(ip) :: nshe
        integer(ip) :: nsho
        integer(ip) :: nte
        integer(ip) :: nto
        real(wp) :: pe
        real(wp) :: po
        real(wp) :: toe
        real(wp) :: tusl
        real(wp) :: ze
        real(wp) :: zo
        real(wp) :: summation, dtheta, v(1, 1), a1, b1, c1
        real(wp) :: cp(idp), work(idp), wx(idp), s(idp+1), &
            e(idp), thet(idp), xx(idp), z(idp), u(idp, idp), &
            we(idp, idp, 2), ped(idp, idp, 2), a(4*idp), b(2*idp), &
            wo(idp, idp, 2), pod(idp, idp, 2)

        dimension pe(idp, idp, 2), po(idp, idp, 2), ze(idp, idp, 2), &
            zo(idp, idp, 2), &
            ipse(idp, 2), jzse(idp, 2), ipso(idp, 2), jzso(idp, 2), &
            nshe(2), nsho(2)

        type(SpherepackUtility) :: util

        ns2 = nlat/2
        modn = nlat-ns2-ns2
        nte = (nlat + 1)/2
        nto = nlat-nte
        tusl = ZERO
        toe = ZERO

        ! Compute grid distribution
        dtheta = pi/(nlat-1)
        do i=1, nte
            thet(i) = real(i-1, kind=wp)*dtheta
        end do

        ! Compute weight matrices for even functions
        do mp1=1, 2
            m = mp1-1
            mrank = nlat-m-m
            nem = (mrank+1)/2
            do j=1, nem
                n = 2*j+m-2
                call util%compute_fourier_coefficients(m, n, cp)
                do i=1, nte
                    call util%compute_legendre_polys_from_fourier_coeff(m, n, thet(i), cp, ped(i, j, mp1))
                end do
                if (m>0) ped(1, j, mp1) = ZERO
            end do

            call singular_value_decomposition(ped(m+1, 1, mp1), idp, nem, nem, s, e, u, &
                idp, v(1, 1), idp, work, 10, info)

            do j=1, nem
                s(j) = ONE/(s(j)**2)
            end do

            ! Compute weight matrix as u  s sup -2 u transpose
            do j=1, nte
                do i=1, nte
                    we(i, j, mp1) = ZERO
                end do
            end do
            do i=1, nem
                do j=1, nem
                    summation = ZERO
                    do k=1, nem
                        summation = summation+s(k)*u(i, k)*u(j, k)
                    end do
                    we(i+m, j+m, mp1) = summation
                end do
            end do
        end do

        we(1, 1, 2) = ONE

        ! Compute n**2 basis (even functions)
        do n=1, nlat+nlat-2
            dfn = n
            a(n) = sqrt(dfn*(dfn+ONE))
        end do
        do n=1, nlat-1
            dfn = n
            b(n) = sqrt((dfn+dfn+3.0)/(dfn+dfn-ONE))
        end do
        !
        mxtr = min(nlat-1, nlon/2, mtrunc)
        iip = 2
        do mp1=1, mxtr+1
            m = mp1-1
            iip = 3-iip
            ms2 = mp1/2
            nrank = ms2+ms2
            mrank = nlat-nrank
            nem = (mrank+1)/2

            ! Compute associated legendre functions
            if (m <= 1) then
                do j=1, nem
                    n = 2*j+m-2
                    call util%compute_fourier_coefficients(m, n, cp)
                    do i=1, nte
                        call util%compute_legendre_polys_from_fourier_coeff(m, n, thet(i), cp, ped(i, j+ms2, iip))
                    end do
                    if (m > 0) ped(1, j+ms2, iip) = ZERO
                end do
            else
                do j=1, nem
                    n = 2*j+m-2
                    if (m > 1 .and. n > mxtr) then
                        do i=1, nte
                            u(i, j+ms2) = ped(i, j+ms2, iip)
                        end do
                    else
                        a1 = b(n-1)*a(n+m-3)/a(n+m-1)
                        b1 = a(n-m+1)/a(n+m-1)
                        if (n-m <= 1) then
                            do i=1, nte
                                u(i, j+ms2) = a1*ped(i, j+ms2-1, iip) &
                                    - b1*ped(i, j+ms2, iip)
                            end do
                        else
                            c1 = b(n-1)*a(n-m-1)/a(n+m-1)
                            do i=1, nte
                                u(i, j+ms2) = a1*ped(i, j+ms2-1, iip) &
                                    - b1*ped(i, j+ms2, iip) + c1*u(i, j+ms2-1)
                            end do
                        end if
                    end if
                end do
                do j=1, nem
                    do i=1, nte
                        ped(i, j+ms2, iip) = u(i, j+ms2)
                    end do
                end do
            end if

            if (.not.(ms2 <= 0. .or. nte <= ms2)) then

                ! initialize array with random numbers
                call random_seed()
                call random_number(xx(1:nte))
                it = 0
                even_iteration: do
                    do i=1, nte
                        z(i) = ZERO
                        wx(i) = ZERO
                        do j=1, nte
                            wx(i) = wx(i)+we(i, j, iip)*xx(j)
                        end do
                    end do

                    do j=1, nte
                        if (j == ms2) cycle
                        call accumulate_inner_products(nte, wx, ped(1, j, iip), z)
                    end do

                    do i=1, nte
                        xx(i) = xx(i)-z(i)
                    end do

                    call compute_normal_regular_grid(nte, xx, idp, we(1, 1, iip))

                    it = it+1
                    if (it > 2) exit even_iteration
                end do even_iteration
                do i=1, nte
                    ped(i, ms2, iip) = xx(i)
                end do
            end if
        end do

        ! Reorder if mtrunc is less than nlat-1
        ! case of even functions
        if (modn == 0) then
            nshe(1) = (nlat-mtrunc-1)/2
            nshe(2) = (nlat-mtrunc-2)/2
        else
            nshe(1) = (nlat-mtrunc)/2
            nshe(2) = (nlat-mtrunc-1)/2
        end if

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

        call truncate(0, nte, idp, ped(1, 1, 1), nte, ipse(1, 1))
        call truncate(0, nte, idp, ped(1, 1, 2), nte, ipse(1, 2))

        ! Compute the analysis matrices
        do iip=1, 2
            do i=1, nte
                lock = 0
                do j=1, nte
                    summation = ZERO
                    do k=1, nte
                        summation = summation+ped(k, i, iip)*we(k, j, iip)
                    end do
                    pe(i, j, iip) = ped(i, j, iip)
                    ze(j, i, iip) =  summation
                    if (abs(summation)>MACHINE_EPSILON .and. lock == 0) then
                        lock = 1
                        jzse(i, iip) = j
                    end if
                end do
            end do
        end do

        ! Compute weight matrices for odd functions
        do mp1=1, 2
            m = mp1-1
            mrank = nlat-m-m
            nem = (mrank+1)/2
            nom = mrank-nem
            do j=1, nom
                n = 2*j+m-1
                call util%compute_fourier_coefficients(m, n, cp)
                do i=1, nte
                    call util%compute_legendre_polys_from_fourier_coeff(m, n, thet(i), cp, pod(i, j, mp1))
                end do
                if (modn == 1) pod(nte, j, mp1) = ZERO
            end do
            call singular_value_decomposition(pod(m+1, 1, mp1), idp, nom, nom, s, e, u, &
                idp, v(1, 1), idp, work, 10, info)

            do j=1, nom
                s(j) = ONE/(s(j)**2)
            end do

            ! Compute weight matrix as u  s sup -2 u transpose
            do j=1, nte
                do i=1, nte
                    wo(i, j, mp1) = ZERO
                end do
            end do
            do i=1, nom
                do j=1, nom
                    summation = ZERO
                    do k=1, nom
                        summation = summation+s(k)*u(i, k)*u(j, k)
                    end do
                    wo(i+m, j+m, mp1) = summation
                end do
            end do
        end do

        wo(1, 1, 2) = ONE
        if (modn == 1) then
            wo(nte, nte, 1) = ONE
            wo(nte, nte, 2) = ONE
        end if

        ! Compute n**2 basis (odd functions)
        iip = 2
        do mp1=1, mxtr+1
            iip = 3-iip
            m = mp1-1
            ms2 = mp1/2
            nrank = ms2+ms2
            mrank = nlat-nrank
            nem = (mrank+1)/2
            nom = mrank-nem

            ! Compute associated legendre functions
            if (m <= 1) then
                do j=1, nom
                    n = 2*j+m-1
                    call util%compute_fourier_coefficients(m, n, cp)
                    do i=1, nte
                        call util%compute_legendre_polys_from_fourier_coeff(m, n, thet(i), cp, pod(i, j+ms2, iip))
                    end do
                    if (modn == 1) pod(nte, j+ms2, iip) = ZERO
                    if (m>0) pod(1, j+ms2, iip) = ZERO
                end do
            else
                do j=1, nom
                    n = 2*j+m-1
                    if (m > 1 .and. n > mxtr) then
                        do i=1, nte
                            u(i, j+ms2) = pod(i, j+ms2, iip)
                        end do
                    else
                        a1 = b(n-1)*a(n+m-3)/a(n+m-1)
                        b1 = a(n-m+1)/a(n+m-1)
                        if (n-m <= 1) then
                            do i=1, nte
                                u(i, j+ms2) = a1*pod(i, j+ms2-1, iip) &
                                    - b1*pod(i, j+ms2, iip)
                            end do
                        else
                            c1 = b(n-1)*a(n-m-1)/a(n+m-1)
                            do i=1, nte
                                u(i, j+ms2) = a1*pod(i, j+ms2-1, iip) &
                                    - b1*pod(i, j+ms2, iip) + c1*u(i, j+ms2-1)
                            end do
                        end if
                    end if
                    if (modn == 1) u(nte, j+ms2) = ZERO
                end do
                do j=1, nom
                    do i=1, nte
                        pod(i, j+ms2, iip) = u(i, j+ms2)
                    end do
                end do
            end if

            if (.not.(ms2 <= 0. .or. nto <= ms2)) then

                ! initialize array with random numbers
                call random_number(xx(1:nte))

                if (modn == 1) xx(nte) = ZERO

                it = 0
                odd_iteration: do
                    do i=1, nte
                        z(i) = ZERO
                        wx(i) = ZERO
                        do j=1, nto
                            wx(i) = wx(i)+wo(i, j, iip)*xx(j)
                        end do
                    end do

                    do j=1, nto
                        if (j == ms2) cycle
                        call accumulate_inner_products(nte, wx, pod(1, j, iip), z(1))
                    end do

                    do i=1, nte
                        xx(i) = xx(i)-z(i)
                    end do

                    call compute_normal_regular_grid(nte, xx, idp, wo(1, 1, iip))

                    it = it+1
                    if (it > 2) exit odd_iteration
                end do odd_iteration

                do i=1, nte
                    pod(i, ms2, iip) = xx(i)
                end do

                if (modn == 1) pod(nte, ms2, iip) = ZERO
            end if
        end do

        ! Reorder if mtrunc is less than nlat-1
        ! case of odd functions
        if (modn == 0) then
            nsho(1) = (nlat-mtrunc)/2
            nsho(2) = (nlat-mtrunc-1)/2
        else
            nsho(1) = (nlat-mtrunc-1)/2
            nsho(2) = (nlat-mtrunc-2)/2
        end if

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

        call truncate(0, nte, idp, pod(1, 1, 1), nto, ipso(1, 1))
        call truncate(0, nte, idp, pod(1, 1, 2), nto, ipso(1, 2))

        ! Compute the analysis matrices (odd functions)
        do iip=1, 2
            do i=1, nto
                lock = 0
                do j=1, nto
                    summation = ZERO
                    do k=1, nte
                        summation = summation+pod(k, i, iip)*wo(k, j, iip)
                    end do
                    po(i, j, iip) = pod(i, j, iip)
                    zo(j, i, iip) = summation
                    if (abs(summation)>MACHINE_EPSILON .and. lock == 0) then
                        lock = 1
                        jzso(i, iip) = j
                    end if
                end do
            end do
        end do

    end subroutine shpei_lower_utility_routine

    subroutine shpe_lower_utility_routine(nlat, nlon, isym, mtrunc, sx, sy, idxy, ierror, &
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
        integer(ip) :: m
        integer(ip) :: modn
        integer(ip) :: mp1
        integer(ip) :: mpm
        integer(ip) :: mrank
        integer(ip) :: ms2
        integer(ip) :: mtrunc
        integer(ip) :: mxtr
        integer(ip) :: nec
        integer(ip) :: nem
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: noc
        integer(ip) :: nom
        integer(ip) :: nrank
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

        dimension sx(idxy, nlon), sy(idxy, nlon), nshe(2), nsho(2), &
            pe(idp, idp, 2), po(idp, idp, 2), ze(idp, idp, 2), zo(idp, idp, 2), &
            ipse(idp, 2), jzse(idp, 2), ipso(idp, 2), jzso(idp, 2), &
            xe(idp, 2), xo(idp, 2), ye(idp, 2), yo(idp, 2)

        ns2 = nlat/2
        modn = nlat-ns2-ns2
        nte = (nlat + 1)/2
        nto = nlat-nte

        if (modn == 0) then
            nshe(1) = (nlat-mtrunc-1)/2
            nshe(2) = (nlat-mtrunc-2)/2
            nsho(1) = (nlat-mtrunc)/2
            nsho(2) = (nlat-mtrunc-1)/2
        else
            nshe(1) = (nlat-mtrunc)/2
            nshe(2) = (nlat-mtrunc-1)/2
            nsho(1) = (nlat-mtrunc-1)/2
            nsho(2) = (nlat-mtrunc-2)/2
        end if
        mxtr = min(nlat-1, nlon/2, mtrunc)
        iip = 2

        main_loop: do mp1=1, mxtr+1
            iip = 3-iip
            if (mxtr == nlat-1 .and. mp1 <= 2) then
                do i=1, nlat
                    sy(i, mp1) = sx(i, mp1)
                end do
                if (mp1 == 2) then
                    sy(1, 2) = ZERO
                    sy(nlat, 2) = ZERO
                end if
                if (3 <= nlon) then
                    sy(1, 3) = ZERO
                    sy(nlat, 3) = ZERO
                    do i=2, nlat-1
                        sy(i, 3) = sx(i, 3)
                    end do
                end if
                cycle main_loop
            end if

            m = mp1-1
            mpm = max(1, m+m)
            ms2 = mp1/2
            mrank = min(nlat-m, nlat-ms2-ms2)
            nrank = nlat-mrank
            nem = (mrank+1)/2-nshe(iip)
            nom = mrank-(mrank+1)/2-nsho(iip)
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
            if (3*nec<2*nem .or. nem == 0) then
                call matrix_multiplication(nte, nec, idp, pe(1, 1, iip), nte, idp, &
                    ze(1, 1, iip), xe, ye, ipse(1, iip), jzse(1, iip))
                do i=1, nte
                    ye(i, 1) = xe(i, 1)-ye(i, 1)
                end do
                if (mpm<nlon .and. m /= 0) then
                    do i=1, nte
                        ye(i, 2) = xe(i, 2)-ye(i, 2)
                    end do
                end if
            else
                call matrix_multiplication(nte, nem, idp, pe(1, nec+1, iip), nte, idp, &
                    ze(1, nec+1, iip), xe, ye, ipse(nec+1, iip), jzse(nec+1, iip))
            end if
            if (3*noc<2*nom .or. nom == 0) then
                call matrix_multiplication(nto, noc, idp, po(1, 1, iip), nto, idp, &
                    zo(1, 1, iip), xo, yo, ipso(1, iip), jzso(1, iip))
                do i=1, nte
                    yo(i, 1) = xo(i, 1)-yo(i, 1)
                end do
                if (mpm<nlon .and. m /= 0) then
                    do i=1, nte
                        yo(i, 2) = xo(i, 2)-yo(i, 2)
                    end do
                end if
            else
                call matrix_multiplication(nto, nom, idp, po(1, noc+1, iip), nto, idp, &
                    zo(1, noc+1, iip), xo, yo, ipso(noc+1, iip), jzso(noc+1, iip))
            end if
            do i=1, nte
                sy(i, mpm) = ye(i, 1)+yo(i, 1)
                sy(nlat+1-i, mpm) = ye(i, 1)-yo(i, 1)
            end do
            if (mpm<nlon .and. m /= 0) then
                do i=1, nte
                    sy(i, mpm+1) = ye(i, 2)+yo(i, 2)
                    sy(nlat+1-i, mpm+1) = ye(i, 2)-yo(i, 2)
                end do
            end if

        end do main_loop

        js = mxtr+mxtr+2

        do j=js, nlon
            do i=1, nlat
                sy(i, j) = ZERO
            end do
        end do

    end subroutine shpe_lower_utility_routine

    subroutine matrix_multiplication(lr, lc, ld, a, mc, md, b, x, y, is, js)

        real(wp) :: a(ld, *), b(md, *)
        integer(ip) :: i
        integer(ip) :: is(*)
        integer(ip) :: j
        integer(ip) :: js(*)
        integer(ip) :: k
        integer(ip) :: kmx
        integer(ip) :: lc
        integer(ip) :: ld
        integer(ip) :: lr
        integer(ip) :: mc
        integer(ip) :: md
        real(wp) :: sum1
        real(wp) :: sum2
        real(wp) :: x(ld, 2), y(ld, 2)

        kmx = min(lr+1, ld)
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

    end subroutine matrix_multiplication

    subroutine compute_normal_regular_grid(n, x, id, q)

        ! Dummy arguments
        integer(ip), intent(in)    :: n
        real(wp),    intent(inout) :: x(n)
        integer(ip), intent(in)    :: id
        real(wp),    intent(in)    :: q(id, n)

        ! Local variables
        integer(ip) :: i, j
        real(wp)    :: summation, sqs

        ! Normalize x
        sqs = ZERO
        do i=1, n
            summation = ZERO
            do j=1, n
                summation = summation + q(i, j)*x(j)
            end do
            sqs = sqs+summation*x(i)
        end do

        x = x/sqrt(sqs)

    end subroutine compute_normal_regular_grid

    ! Purpose:
    !
    ! Performs the singular value decomposition of a real rectangular matrix.
    !
    !    This routine reduces an m by n matrix a to diagonal form by orthogonal
    !    transformations u and v.  The diagonal elements s(i) are the singular
    !    values of a.  The columns of u are the corresponding left singular
    !    vectors, and the columns of v the right singular vectors.
    !
    !    The form of the singular value decomposition is then
    !
    !      a(mxn) = u(mxm) * s(mxn) * transpose(v(nxn))
    !
    !  Reference:
    !
    !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart, 
    !    LINPACK User's Guide, 
    !    SIAM, 1979, 
    !    ISBN13: 978-0-898711-72-1, 
    !    LC: QA214.L56.
    !
    !  Parameters:
    !
    !    Input/output, real(wp) a(lda, n).  on input, the m by n
    !    matrix whose singular value decomposition is to be computed.
    !    on output, the matrix has been destroyed.  depending on the user's
    !    requests, the matrix may contain other useful information.
    !
    !    input, integer(ip) lda, the leading dimension of the array a.
    !    lda must be at least n.
    !
    !    input, integer(ip) m, the number of rows of the matrix.
    !
    !    input, integer(ip) n, the number of columns of the matrix a.
    !
    !    output, real(wp) s(mm), where mm = max(m+1, n).  the first
    !    min(m, n) entries of s contain the singular values of a arranged in
    !    descending order of magnitude.
    !
    !    output, real(wp) e(mm), where mm = max(m+1, n).  ordinarily
    !    contains zeros.  however see the discussion of info for exceptions.
    !
    !    output, real(wp) u(ldu, k).  if joba = 1 then k = m;
    !    if 2 <= joba, then k = min(m, n).  u contains the m by m matrix of
    !    left singular vectors.  u is not referenced if joba = 0.  if m <= n
    !    or if joba = 2, then u may be identified with a in the subroutine call.
    !
    !    input, integer(ip) ldu, the leading dimension of the array u.
    !    ldu must be at least m.
    !
    !    output, real(wp) v(ldv, n), the n by n matrix of right singular
    !    vectors.  v is not referenced if job is 0.  if n <= m, then v may be
    !    identified with a in the subroutine call.
    !
    !    input, integer(ip) ldv, the leading dimension of the array v.
    !    ldv must be at least n.
    !
    !    workspace, real(wp) work(m).
    !
    !    input, integer(ip) job, controls the computation of the singular
    !    vectors.  it has the decimal expansion ab with the following meaning:
    !      a =  0, do not compute the left singular vectors.
    !      a =  1, return the m left singular vectors in u.
    !      a >= 2, return the first min(m, n) singular vectors in u.
    !      b =  0, do not compute the right singular vectors.
    !      b =  1, return the right singular vectors in v.
    !
    !    output, integer(ip) info, status indicator.
    !    the singular values (and their corresponding singular vectors)
    !    s(info+1), s(info+2), ..., s(mn) are correct.  here mn = min(m, n).
    !    thus if info is 0, all the singular values and their vectors are
    !    correct.  in any event, the matrix b = u' * a * v is the bidiagonal
    !    matrix with the elements of s on its diagonal and the elements of e on
    !    its superdiagonal.  thus the singular values of a and b are the same.
    !
    subroutine singular_value_decomposition(a, lda, m, n, s, e, u, ldu, v, ldv, work, job, info)

        ! Dummy arguments
        real(wp),    intent(inout) :: a(lda, n)
        integer(ip), intent(in)    :: lda
        integer(ip), intent(in)    :: m
        integer(ip), intent(in)    :: n
        real(wp),    intent(out)   :: s(*)
        real(wp),    intent(out)   :: e(*)
        real(wp),    intent(out)   :: u(ldu, m)
        integer(ip), intent(in)    :: ldu
        real(wp),    intent(out)   :: v(ldv, n)
        integer(ip), intent(in)    :: ldv
        real(wp),    intent(out)   :: work(m)
        integer(ip), intent(in)    :: job
        integer(ip), intent(out)   :: info

        ! Local variables
        real(wp) b
        real(wp) c
        real(wp) cs
        real(wp) el
        real(wp) emm1
        real(wp) f
        real(wp) g
        integer(ip) iter
        integer(ip) j
        integer(ip) job_u
        integer(ip) k
        integer(ip) k_case
        integer(ip) kk
        integer(ip) l
        integer(ip) ll
        integer(ip) lls
        integer(ip) ls
        integer(ip) lu
        integer(ip), parameter :: MAX_ITER = 30
        integer(ip) mm
        integer(ip) mm1
        integer(ip) mn
        integer(ip) nct
        integer(ip) nctp1
        integer(ip) ncu
        integer(ip) nrt
        integer(ip) nrtp1
        real(wp) scale_factor
        real(wp) shift
        real(wp) sl
        real(wp) sm
        real(wp) smm1
        real(wp) sn
        real(wp) t
        real(wp) t1
        real(wp) test
        real(wp) ztest
        logical :: u_desired, v_desired

        !  Determine what is to be computed.
        u_desired = .false.
        v_desired = .false.
        job_u = mod(job, 100)/10

        if (1 < job_u) then
            ncu = min(m, n)
        else
            ncu = m
        end if

        if (job_u /= 0) u_desired = .true.

        if (mod(job, 10) /= 0) v_desired = .true.

        !  Reduce a to bidiagonal form, storing the diagonal elements
        !  in s and the super-diagonal elements in e.
        info = 0
        nct = min(m - 1, n)
        nrt = max(0, min(m, n - 2))
        lu = max(nct, nrt)

        do l = 1, lu

            !  Compute the transformation for the l-th column and
            !  place the l-th diagonal in s(l).
            if (l <= nct) then
                s(l) = get_norm2(m-l+1, a(l, l), 1)
                if (s(l) /= ZERO) then
                    if (a(l, l) /= ZERO) s(l) = sign(s(l), a(l, l))
                    call scale_vector_by_constant( m-l+1, ONE / s(l), a(l, l), 1)
                    a(l, l) = ONE + a(l, l)
                end if
                s(l) = -s(l)
            end if

            do j = l + 1, n

                !  Apply the transformation.
                if (l <= nct .and. s(l) /= ZERO) then
                    t = -get_dot_product( m-l+1, a(l, l), 1, a(l, j), 1) / a(l, l)
                    call daxpy(m-l+1, t, a(l, l), 1, a(l, j), 1)
                end if

                !  Place the l-th row of a into e for the
                !  subsequent calculation of the row transformation.
                e(j) = a(l, j)
            end do

            !  Place the transformation in u for subsequent back multiplication.
            if (u_desired .and. l <= nct) u(l:m, l) = a(l:m, l)

            !  Compute the l-th row transformation and place the
            !  l-th superdiagonal in e(l).
            !
            if (l <= nrt) then
                e(l) = get_norm2(n - l, e(l+1), 1)
                if (e(l) /= ZERO) then
                    if (e(l+1) /= ZERO) then
                        e(l) = sign(e(l), e(l+1))
                    end if
                    call scale_vector_by_constant( n-l, ONE / e(l), e(l+1), 1)
                    e(l+1) = ONE + e(l+1)
                end if

                e(l) = -e(l)
                !
                !  Apply the transformation.
                !
                if (l + 1 <= m .and. e(l) /= ZERO) then

                    work(l+1:m) = ZERO

                    do j = l + 1, n
                        call daxpy(m-l, e(j), a(l+1, j), 1, work(l+1), 1)
                    end do

                    do j = l + 1, n
                        call daxpy(m-l, -e(j)/e(l+1), work(l+1), 1, a(l+1, j), 1)
                    end do

                end if

                !  Place the transformation in v for subsequent back multiplication.
                if (v_desired) v(l+1:n, l) = e(l+1:n)
            end if
        end do

        ! Set up the final bidiagonal matrix of order mn.
        mn = min(m + 1, n)
        nctp1 = nct + 1
        nrtp1 = nrt + 1

        if (nct < n) s(nctp1) = a(nctp1, nctp1)

        if (m < mn) s(mn) = ZERO

        if (nrtp1 < mn) e(nrtp1) = a(nrtp1, mn)

        e(mn) = ZERO

        !  If required, generate u.
        if (u_desired) then

            u(1:m, nctp1:ncu) = ZERO

            do j = nctp1, ncu
                u(j, j) = ONE
            end do

            do ll = 1, nct

                l = nct - ll + 1

                if (s(l) /= ZERO) then

                    do j = l + 1, ncu
                        t = -get_dot_product( m-l+1, u(l, l), 1, u(l, j), 1) / u(l, l)
                        call daxpy(m-l+1, t, u(l, l), 1, u(l, j), 1)
                    end do

                    u(l:m, l) = -u(l:m, l)
                    u(l, l) = ONE + u(l, l)
                    u(1:l-1, l) = ZERO

                else
                    u(1:m, l) = ZERO
                    u(l, l) = ONE
                end if
            end do
        end if

        !  If it is required, generate v.
        if (v_desired) then

            do ll = 1, n
                l = n - ll + 1
                if (l <= nrt .and. e(l) /= ZERO) then
                    do j = l + 1, n
                        t = -get_dot_product( n-l, v(l+1, l), 1, v(l+1, j), 1) / v(l+1, l)
                        call daxpy(n-l, t, v(l+1, l), 1, v(l+1, j), 1)
                    end do
                end if
                v(1:n, l) = ZERO
                v(l, l) = ONE
            end do
        end if

        !  Main iteration loop for the singular values.
        mm = mn
        iter = 0
        main_iteration_loop: do

            if (mn <= 0) exit main_iteration_loop

            !  If too many iterations have been performed, set flag and return.

            if (MAX_ITER <= iter) then
                info = mn
                return
            end if
            !
            !  This section of the program inspects for
            !  negligible elements in the s and e arrays.
            !
            !  On completion the variables k_case and l are set as follows:
            !
            !  k_case = 1     if s(mn) and e(l-1) are negligible and l < mn
            !  k_case = 2     if s(l) is negligible and l < mn
            !  k_case = 3     if e(l-1) is negligible, l < mn, and
            !               s(l), ..., s(mn) are not negligible (qr step).
            !  k_case = 4     if e(mn-1) is negligible (convergence).
            !
            do ll = 1, mn
                l = mn - ll
                if (l == 0) exit
                test = abs(s(l)) + abs(s(l+1))
                ztest = test + abs(e(l))
                if (ztest == test) then
                    e(l) = ZERO
                    exit
                end if
            end do

            if (l == mn - 1) then
                k_case = 4
            else
                do lls = l + 1, mn + 1
                    ls = mn - lls + l + 1

                    if (ls == l) exit

                    test = ZERO

                    if (ls /= mn) test = test + abs(e(ls))

                    if (ls /= l + 1) test = test + abs(e(ls-1))

                    ztest = test + abs(s(ls))

                    if (ztest == test) then
                        s(ls) = ZERO
                        exit
                    end if
                end do

                if (ls == l) then
                    k_case = 3
                else if (ls == mn) then
                    k_case = 1
                else
                    k_case = 2
                    l = ls
                end if
            end if

            l = l + 1

            !  Deflate negligible s(mn).
            select case (k_case)
                case (1)
            		
                    mm1 = mn - 1
                    f = e(mn-1)
                    e(mn-1) = ZERO
            		
                    do kk = l, mm1
            		
                        k = mm1 - kk + l
                        t1 = s(k)
                        call construct_givens_plane_rotation(t1, f, cs, sn)
                        s(k) = t1
            		
                        if (k /= l) then
                            f = -sn * e(k-1)
                            e(k-1) = cs * e(k-1)
                        end if
            		
                        if (v_desired) call apply_plane_rotation(n, v(1, k), 1, v(1, mn), 1, cs, sn)
                    end do
                case (2)
            		
                    f = e(l-1)
                    e(l-1) = ZERO
            		
                    do k = l, mn
                        t1 = s(k)
                        call construct_givens_plane_rotation(t1, f, cs, sn)
                        s(k) = t1
                        f = -sn * e(k)
                        e(k) = cs * e(k)
                        if (u_desired) call apply_plane_rotation(m, u(1, k), 1, u(1, l-1), 1, cs, sn)
                    end do
                case (3)

                    !  Calculate the shift.
                    scale_factor = max(abs(s(mn)), abs(s(mn - 1)), abs(e(mn - 1)), &
                        abs(s(l)), abs(e(l)))
            		
                    sm = s(mn)/scale_factor
                    smm1 = s(mn-1)/scale_factor
                    emm1 = e(mn-1)/scale_factor
                    sl = s(l)/scale_factor
                    el = e(l)/scale_factor
                    b = ((smm1 + sm) * (smm1 - sm) + emm1**2)/TWO
                    c = (sm**2) * (emm1**2)
                    shift = ZERO
            		
                    if (b /= ZERO .or. c /= ZERO) then
                        shift = sqrt(b**2 + c)
                        if (b < ZERO) shift = -shift
                        shift = c / (b + shift)
                    end if
            		
                    f = ( sl + sm) * ( sl - sm) + shift
                    g = sl * el

                    !  Chase zeros.
                    mm1 = mn - 1
            		
                    do k = l, mm1
            		
                        call construct_givens_plane_rotation(f, g, cs, sn)
            		
                        if (k /= l) e(k-1) = f
            		
                        f = cs * s(k) + sn * e(k)
                        e(k) = cs * e(k) - sn * s(k)
                        g = sn * s(k+1)
                        s(k+1) = cs * s(k+1)
            		
                        if (v_desired) call apply_plane_rotation(n, v(1, k), 1, v(1, k+1), 1, cs, sn)

                        call construct_givens_plane_rotation(f, g, cs, sn)

                        s(k) = f
                        f = cs * e(k) + sn * s(k+1)
                        s(k+1) = -sn * e(k) + cs * s(k+1)
                        g = sn * e(k+1)
                        e(k+1) = cs * e(k+1)
            		
                        if (u_desired .and. k < m) call apply_plane_rotation(m, u(1, k), 1, u(1, k+1), 1, cs, sn)

                    end do
            		
                    e(mn-1) = f
                    iter = iter + 1
                case (4)

                    !  Make the singular value nonnegative.
                    if (s(l) < ZERO) then
                        s(l) = -s(l)
                        if (v_desired) v(1:n, l) = -v(1:n, l)
                    end if

                    !  Order the singular value.
                    do
            		
                        if (l == mm) exit
            		
                        if (s(l+1) <= s(l)) exit
            		
                        t = s(l)
                        s(l) = s(l+1)
                        s(l+1) = t
            		
                        if (v_desired .and. l < n) call swap_vectors( n, v(1, l), 1, v(1, l+1), 1)

                        if (u_desired .and. l < m) call swap_vectors(m, u(1, l), 1, u(1, l+1), 1)
                        l = l + 1
                    end do
                    iter = 0
                    mn = mn - 1
            end select
        end do main_iteration_loop

    end subroutine singular_value_decomposition

    ! Purpose:
    !
    ! Computes constant times a vector plus a vector.
    !
    ! Jack dongarra, linpack, 3/11/78.
    ! Modified 12/3/93, array(1) declarations changed to array(*)
    ! Modified 01/27/17, whole array operations to aid compiler optimization
    !
    subroutine daxpy(n, da, dx, incx, dy, incy)

        ! Dummy arguments
        integer(ip), intent(in)    :: n
        real(wp),    intent(in)    :: da
        real(wp),    intent(in)    :: dx(*)
        integer(ip), intent(in)    :: incx
        real(wp),    intent(inout) :: dy(*)
        integer(ip), intent(in)    :: incy

        ! Local variables
        integer(ip) :: i, ix, iy, m

        if (n <= 0) then
            return
        else if (da == ZERO) then
            return
        else if (incx /= 1 .or. incy /= 1) then

            !  Code for unequal increments or equal increments
            !  not equal to 1.
            if (0 <= incx) then
                ix = 1
            else
                ix = (-n + 1) * incx + 1
            end if

            if (0 <= incy) then
                iy = 1
            else
                iy = (-n + 1) * incy + 1
            end if

            do i = 1, n
                dy(iy) = dy(iy) + da * dx(ix)
                ix = ix + incx
                iy = iy + incy
            end do
        else
            !  Code for both increments equal to 1.
            m = mod(n, 4)

            dy(1:m) = dy(1:m) + da * dx(1:m)

            do i = m + 1, n, 4
                dy(i) = dy(i) + da * dx(i)
                dy(i+1) = dy(i+1) + da * dx(i+1)
                dy(i+2) = dy(i+2) + da * dx(i+2)
                dy(i+3) = dy(i+3) + da * dx(i+3)
            end do
        end if

    end subroutine daxpy

    ! Purpose:
    !
    ! Forms the dot product of two vectors.
    ! Jack dongarra, linpack, 3/11/78.
    ! Modified 12/3/93, array(1) declarations changed to array(*)
    ! Modified 01/27/17, the function now wraps around the intrinsic dot_product
    !
    pure function get_dot_product(n, dx, incx, dy, incy) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: n
        real(wp),    intent(in) :: dx(*)
        integer(ip), intent(in) :: incx
        real(wp),    intent(in) :: dy(*)
        integer(ip), intent(in) :: incy
        real(wp)                :: return_value

        ! Local variables
        integer(ip) :: x1, xn, xi
        integer(ip) :: y1, yn, yi

        if (0 < incx) then
            x1 = 1
            xn = 1 + (n - 1) * incx
            xi = incx
        else
            x1 = 1 + (n - 1) * incx
            xn = 1
            xi = - incx
        end if

        if (0 < incy) then
            y1 = 1
            yn = 1 + (n - 1) * incy
            yi = incy
        else
            y1 = 1 + (n - 1) * incy
            yn = 1
            yi = - incy
        end if

        !  Let the intrinsic function dot_product take care of optimization.
        return_value = dot_product(dx(x1:xn:xi), dy(y1:yn:yi))

    end function get_dot_product

    ! Purpose:
    !
    !  Returns the euclidean norm of a vector via the function
    !  name, so that
    !
    !     return_value := sqrt( transpose(x) * x)
    !
    ! This version written on 25-October-1982.
    ! Modified on 14-October-1993 to inline the call to DLASSQ.
    ! Sven Hammarling, Nag Ltd.
    !
    ! Modified 01/27/17, function now wraps around the intrinsic norm2
    !
    pure function get_norm2(n, x, incx) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: n
        real(wp),    intent(in) :: x(*)
        integer(ip), intent(in) :: incx
        real(wp)                :: return_value

        ! Local variables
        integer(ip) :: x1, xn, xi

        if (0 < incx) then
            x1 = 1
            xn = 1 + (n - 1) * incx
            xi = incx
        else
            x1 = 1 + (n - 1) * incx
            xn = 1
            xi = -incx
        end if

        !  Let the intrinsic function norm2 take care of optimization.
        return_value = norm2(x(x1:xn:xi))

    end function get_norm2

    ! Purpose:
    !
    ! Applies a plane rotation.
    ! Jack dongarra, linpack, 3/11/78.
    ! Modified 12/3/93, array(1) declarations changed to array(*)
    !
    subroutine apply_plane_rotation(n, x, incx, y, incy, c, s)

        ! Dummy arguments
        integer(ip), intent(in)    :: n
        real(wp),    intent(inout) :: x(*)
        integer(ip), intent(in)    :: incx
        real(wp),    intent(inout) :: y(*)
        integer(ip), intent(in)    :: incy
        real(wp),    intent(in)    :: c
        real(wp),    intent(in)    :: s

        ! Local variables
        integer(ip) :: i, ix, iy
        real(wp)    :: temp

        if (n <= 0) then
            return
        else if (incx == 1 .and. incy == 1) then
            do i = 1, n
                temp = c * x(i) + s * y(i)
                y(i) = c * y(i) - s * x(i)
                x(i) = temp
            end do
        else
            if (0 <= incx) then
                ix = 1
            else
                ix = (-n + 1) * incx + 1
            end if

            if (0 <= incy) then
                iy = 1
            else
                iy = (-n + 1) * incy + 1
            end if

            do i = 1, n
                temp = c * x(ix) + s * y(iy)
                y(iy) = c * y(iy) - s * x(ix)
                x(ix) = temp
                ix = ix + incx
                iy = iy + incy
            end do
        end if

    end subroutine apply_plane_rotation

    ! Purpose:
    !
    ! Constructs a Givens plane rotation.
    !
    ! Given values a and b, this routine computes
    !
    !    sigma = sign(a) if abs(a) >  abs(b)
    !          = sign(b) if abs(a) <= abs(b);
    !
    !    r     = sigma * (a**2 + b**2);
    !
    !    c = a / r if r is not 0
    !      = 1     if r is 0;
    !
    !    s = b / r if r is not 0, 
    !        0     if r is 0.
    !
    !    The computed numbers then satisfy the equation
    !
    !    (  c  s) ( a) = ( r)
    !    ( -s  c) ( b) = ( 0)
    !
    !    The routine also computes
    !
    !    z = s     if abs(a) > abs(b), 
    !      = 1 / c if abs(a) <= abs(b) and c is not 0, 
    !      = 1     if c is 0.
    !
    !    The single value z encodes c and s, and hence the rotation:
    !
    !    if z = 1, set c = 0 and s = 1;
    !    if abs(z) < 1, set c = sqrt(1 - z**2) and s = z;
    !    if abs(z) > 1, set c = 1/ z and s = sqrt(1 - c**2);
    !
    !  Reference:
    !
    !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart, 
    !    LINPACK User's Guide, 
    !    SIAM, 1979, 
    !    ISBN13: 978-0-898711-72-1, 
    !    LC: QA214.L56.
    !
    !    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh, 
    !    Algorithm 539, 
    !    Basic Linear Algebra Subprograms for Fortran Usage, 
    !    ACM Transactions on Mathematical Software, 
    !    Volume 5, Number 3, September 1979, pages 308-323.
    !
    !  Dummy arguments:
    !
    !    input/output, real sa, sb.  on input, sa and sb are the values
    !    a and b.  on output, sa is overwritten with r, and sb is
    !    overwritten with z.
    !
    !    Output, real c, s, the cosine and sine of the
    !    Givens rotation.
    !
    subroutine construct_givens_plane_rotation(sa, sb, c, s)

        ! Dummy arguments
        real(wp), intent(inout) :: sa
        real(wp), intent(inout) :: sb
        real(wp), intent(out) :: s
        real(wp), intent(out) :: c

        ! Local variables
        real(wp) :: scale_factor, r, z, roe

        if (abs(sb) < abs(sa)) then
            roe = sa
        else
            roe = sb
        end if

        scale_factor = abs(sa) + abs(sb)

        if (scale_factor == ZERO) then
            c = ONE
            s = ZERO
            r = ZERO
        else
            r = scale_factor * hypot((sa/scale_factor), (sb/scale_factor))
            r = sign(ONE, roe) * r
            c = sa/r
            s = sb/r
        end if

        if (ZERO < abs(c) .and. abs(c) <= s) then
            z = ONE/c
        else
            z = s
        end if

        sa = r
        sb = z

    end subroutine construct_givens_plane_rotation

    ! Purpose:
    !
    ! Scales a vector by a constant.
    !
    ! Jack dongarra, linpack, 3/11/78.
    ! Modified 3/93 to return if incx <= 0.0
    ! Modified 12/3/93, array(1) declarations changed to array(*)
    ! Modified 01/27/17, uses array operations to aid compiler optimization
    !
    subroutine scale_vector_by_constant(n, sa, x, incx)

        ! Dummy arguments
        integer(ip), intent(in)    :: n
        real(wp),    intent(in)    :: sa
        real(wp),    intent(inout) :: x(*)
        integer(ip), intent(in)    :: incx

        ! Local variables
        integer(ip) :: i, ix, m

        if (n <= 0) then
            return
        else if (incx == 1) then

            m = mod(n, 5)

            x(1:m) = sa * x(1:m)

            do i = m+1, n, 5
                x(i) = sa * x(i)
                x(i+1) = sa * x(i+1)
                x(i+2) = sa * x(i+2)
                x(i+3) = sa * x(i+3)
                x(i+4) = sa * x(i+4)
            end do
        else
            if (0 <= incx) then
                ix = 1
            else
                ix = (-n + 1) * incx + 1
            end if

            do i = 1, n
                x(ix) = sa * x(ix)
                ix = ix + incx
            end do

        end if

    end subroutine scale_vector_by_constant

    ! Purpose:
    !
    ! Interchanges two vectors.
    ! uses unrolled loops for increments equal one.
    ! Jack dongarra, linpack, 3/11/78.
    ! Modified 12/3/93, array(1) declarations changed to array(*)
    ! Modified 01/27/17, whole array operations to improve compiler optimization
    !
    subroutine swap_vectors(n, x, incx, y, incy)

        ! Dummy arguments
        integer(ip), intent(in)    :: n
        real(wp),    intent(inout) :: x(*)
        integer(ip), intent(in)    :: incx
        real(wp),    intent(inout) :: y(*)
        integer(ip), intent(in)    :: incy

        ! Local variables
        integer(ip) :: m, i, ix, iy
        real(wp)    :: temp

        if (n <= 0) then
            return
        else if (incx == 1 .and. incy == 1) then
            m = mod(n, 3)
            do i = 1, m
                temp = x(i)
                x(i) = y(i)
                y(i) = temp
            end do

            do i = m + 1, n, 3
                temp = x(i)
                x(i) = y(i)
                y(i) = temp

                temp = x(i+1)
                x(i+1) = y(i+1)
                y(i+1) = temp

                temp = x(i+2)
                x(i+2) = y(i+2)
                y(i+2) = temp
            end do
        else

            if (0 <= incx) then
                ix = 1
            else
                ix = (-n + 1) * incx + 1
            end if

            if (0 <= incy) then
                iy = 1
            else
                iy = (-n + 1) * incy + 1
            end if

            do i = 1, n
                temp = x(ix)
                x(ix) = y(iy)
                y(iy) = temp
                ix = ix + incx
                iy = iy + incy
            end do
        end if

    end subroutine swap_vectors

end submodule scalar_projection_regular_grid
