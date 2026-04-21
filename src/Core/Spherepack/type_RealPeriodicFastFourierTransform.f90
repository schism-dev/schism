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
!     This file contains a multiple fft package for spherepack. It
!     includes code and documentation for performing fast Fourier
!     transforms (see subroutines hrffti, hrfftf and hrfftb)
!
module type_RealPeriodicFastFourierTransform

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI, &
        TWO_PI, &
        odd

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: hrffti, hrfftf, hrfftb
    
    ! Parameters confined to the module
    real(wp),    parameter :: ZERO = 0.0_wp
    real(wp),    parameter :: ONE = 1.0_wp
    real(wp),    parameter :: TWO = 2.0_wp
    real(wp),    parameter :: SQRT_2 = sqrt(TWO)
    integer(ip), parameter :: NUMBER_OF_FACTORS = 15_ip

    type, public :: RealPeriodicFastFourierTransform
    contains
        ! Type-bound procedures
        procedure, nopass :: initialize => hrffti
        procedure, nopass :: forward => hrfftf
        procedure, nopass :: backward => hrfftb
    end type RealPeriodicFastFourierTransform

contains

    ! Purpose:
    !
    !     subroutine hrffti(n, wsave)
    !
    !     subroutine hrffti initializes the array wsave which is used in
    !     both hrfftf and hrfftb. the prime factorization of n together
    !     with a tabulation of the trigonometric functions are computed and
    !     stored in wsave.
    !
    !     input parameter
    !
    !     n       the length of the sequence to be transformed.
    !
    !     output parameter
    !
    !     wsave   a work array which must be dimensioned at least n+15.
    !             the same work array can be used for both hrfftf and
    !             hrfftb as long as n remains unchanged. different wsave
    !             arrays are required for different values of n. the
    !             contents of wsave must not be changed between calls
    !             of hrfftf or hrfftb.
    !
    pure subroutine hrffti(n, wsave)

        ! Dummy arguments
        integer(ip), intent(in)  :: n
        real(wp),    intent(out) :: wsave(n+NUMBER_OF_FACTORS)

        ! Local variables
        integer(ip) :: iw1, iw2

        if (n > 1) then

            ! Set wavetable index pointers
            iw1 = 1
            iw2 = n + 1

            call precompute_factorization_and_trig_lookup_table(n, wsave(iw1:), wsave(iw2:))
        end if

    end subroutine hrffti

    ! Purpose:
    !
    !     subroutine hrfftf(m, n, r, mdimr, wsave)
    !
    !     Computes the Fourier coefficients of m real
    !     perodic sequences (Fourier analysis); i.e. hrfftf computes the
    !     real fft of m sequences each with length n. the transform is
    !     defined below at output parameter r.
    !
    !     input parameters
    !
    !     m       the number of sequences.
    !
    !     n       the length of all m sequences.  the method is most
    !             efficient when n is a product of small primes. n may
    !             change as long as different work arrays are provided
    !
    !     r       r(mdimr, n) is a two dimensional real array that contains mdimr
    !             sequences each with length n.
    !
    !     mdimr   the first dimension of the r array as it appears
    !             in the program that calls hrfftf. mdimr must be
    !             greater than or equal to m.
    !
    !
    !     wsave   a work array with at least least n+15 locations
    !             in the program that calls hrfftf. the wsave array must be
    !             initialized by calling subroutine hrffti(n, wsave) and a
    !             different wsave array must be used for each different
    !             value of n. this initialization does not have to be
    !             repeated so long as n remains unchanged thus subsequent
    !             transforms can be obtained faster than the first.
    !             the same wsave array can be used by hrfftf and hrfftb.
    !
    !     work    a real work array with m*n locations.
    !
    !
    !     output parameters
    !
    !     r      for all j=1, ..., m
    !
    !             r(j, 1) = the sum from i=1 to i=n of r(j, i)
    !
    !             if n is even set l =n/2   , if n is odd set l = (n+1)/2
    !
    !               then for k = 2, ..., l
    !
    !                  r(j, 2*k-2) = the sum from i = 1 to i = n of
    !
    !                       r(j, i)*cos((k-1)*(i-1)*2*pi/n)
    !
    !                  r(j, 2*k-1) = the sum from i = 1 to i = n of
    !
    !                      -r(j, i)*sin((k-1)*(i-1)*2*pi/n)
    !
    !             if n is even
    !
    !                  r(j, n) = the sum from i = 1 to i = n of
    !
    !                       (-1)**(i-1)*r(j, i)
    !
    !      *****  note
    !                  this transform is unnormalized since a call of hrfftf
    !                  followed by a call of hrfftb will multiply the input
    !                  sequence by n.
    !
    !     wsave   contains results which must not be destroyed between
    !             calls of hrfftf or hrfftb.
    !
    subroutine hrfftf(m, n, r, mdimr, wsave)

        ! Dummy arguments
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: n
        real(wp),    intent(inout)  :: r(mdimr, n)
        integer(ip), intent(in)     :: mdimr
        real(wp),    intent(in)     :: wsave(n+NUMBER_OF_FACTORS)

        ! Local variables
        integer(ip) :: iw1, iw2

        if (n > 1) then

            ! Set workspace index pointers
            iw1 = 1
            iw2 = n + 1

            call forward_lower_utility_routine(m, n, r, mdimr, &
                wsave(iw1:), wsave(iw2:))
        end if

    end subroutine hrfftf

    ! Purpose:
    !
    !     subroutine hrfftb(m, n, r, mdimr, wsave, work)
    !
    !     subroutine hrfftb computes the real perodic sequence of m
    !     sequences from their Fourier coefficients (Fourier synthesis).
    !     the transform is defined below at output parameter r.
    !
    !     input parameters
    !
    !     m       the number of sequences.
    !
    !     n       the length of all m sequences.  the method is most
    !             efficient when n is a product of small primes. n may
    !             change as long as different work arrays are provided
    !
    !     r       r(mdimr, n) is a two dimensional real array that contains
    !             the Fourier coefficients of m sequences each with
    !             length n.
    !
    !     mdimr   the first dimension of the r array as it appears
    !             in the program that calls hrfftb. mdimr must be
    !             greater than or equal to m.
    !
    !     wsave   a work array which must be dimensioned at least n+15.
    !             in the program that calls hrfftb. the wsave array must be
    !             initialized by calling subroutine hrffti(n, wsave) and a
    !             different wsave array must be used for each different
    !             value of n. this initialization does not have to be
    !             repeated so long as n remains unchanged thus subsequent
    !             transforms can be obtained faster than the first.
    !             the same wsave array can be used by hrfftf and hrfftb.
    !
    !     work    a real work array with m*n locations.
    !
    !
    !     output parameters
    !
    !     r      for all j=1, ..., m
    !
    !             for n even and for i = 1, ..., n
    !
    !                  r(j, i) = r(j, 1)+(-1)**(i-1)*r(j, n)
    !
    !                       plus the sum from k=2 to k=n/2 of
    !
    !                        2.0*r(j, 2*k-2)*cos((k-1)*(i-1)*2*pi/n)
    !
    !                       -2.0*r(j, 2*k-1)*sin((k-1)*(i-1)*2*pi/n)
    !
    !             for n odd and for i = 1, ..., n
    !
    !                  r(j, i) = r(j, 1) plus the sum from k=2 to k=(n+1)/2 of
    !
    !                       2.0*r(j, 2*k-2)*cos((k-1)*(i-1)*2*pi/n)
    !
    !                      -2.0*r(j, 2*k-1)*sin((k-1)*(i-1)*2*pi/n)
    !
    !      *****  note
    !                  this transform is unnormalized since a call of hrfftf
    !                  followed by a call of hrfftb will multiply the input
    !                  sequence by n.
    !
    !     wsave   contains results which must not be destroyed between
    !             calls of hrfftb or hrfftf.
    !
    subroutine hrfftb(m, n, r, mdimr, wsave)

        ! Dummy arguments
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: n
        real(wp),    intent(inout)  :: r(mdimr, n)
        integer(ip), intent(in)     :: mdimr
        real(wp),    intent(in)     :: wsave(n+NUMBER_OF_FACTORS)

        ! Local variables
        integer(ip) :: iw1, iw2

        if (n > 1) then

            ! Set wavetable index pointers
            iw1 = 1
            iw2 = n + 1

            call backward_lower_utility_routine(m, n, r, mdimr, &
                wsave(iw1:), wsave(iw2:))
        end if

    end subroutine hrfftb

    ! Purpose:
    !
    ! Factors of an integer for floating point computations.
    !
    pure subroutine compute_factorization(n, subtransforms, number_of_factors, integer_factors)

        ! Dummy arguments
        integer(ip), intent(in)  :: n
        integer(ip), intent(in)  :: subtransforms(:)
        integer(ip), intent(out) :: number_of_factors
        real(wp),    intent(out) :: integer_factors(:)

        ! Local variables
        integer(ip) :: nf, ntest, j, factor

        ! Initialize
        factor = 0
        ntest = n
        nf = 0
        j = 0

        select case(n)
            case(:0)
                error stop 'Length n must be a positive integer'
            case(1)
                integer_factors(1) = 1
                number_of_factors = 1
            case default

                block
                    integer(ip) :: factors(size(integer_factors))

                    do while (1 < ntest)

                        ! Increment j
                        j = j + 1

                        ! Choose factor
                        select case (j)
                            case(1:4) ! case(1:size(subtransforms))
                                factor = subtransforms(j)
                            case default
                                factor = factor + 2
                        end select

                        do while ((ntest - factor * (ntest/factor)) == 0)
                            nf = nf + 1
                            factors(nf + 2) = factor
                            ntest = ntest/factor
                            if (factor == 2 .and. nf /= 1) then
                                factors(nf+2:4:(-1)) = factors(nf+1:3:(-1))
                                factors(3) = 2
                            end if
                        end do
                    end do

                    factors(1) = n
                    factors(2) = nf
                    number_of_factors = nf
                    integer_factors = real(factors, kind=wp)
                end block
        end select

    end subroutine compute_factorization

    pure subroutine precompute_factorization_and_trig_lookup_table(n, trig_lookup_table, factors)

        ! Dummy arguments
        integer(ip), intent(in)  :: n
        real(wp),    intent(out) :: trig_lookup_table(:)
        real(wp),    intent(out) :: factors(:)

        ! Local variables
        integer(ip)            :: i, ido, ii, iip, ipm, iis
        integer(ip)            :: j, k1, l1, l2, ld
        integer(ip)            :: nf, nfm1
        integer(ip), parameter :: SUBTRANSFORMS(*) = [4, 2, 3, 5]
        real(wp)               :: theta, d_theta, mesh, product_1

        call compute_factorization(n, SUBTRANSFORMS, nf, factors)

        d_theta = TWO_PI/n
        iis = 0
        nfm1 = nf-1
        l1 = 1

        if (nfm1 /= 0) then
            do k1=1, nfm1
                iip = int(factors(k1+2), kind=ip)
                ld = 0
                l2 = l1*iip
                ido = n/l2
                ipm = iip-1
                do j=1, ipm
                    ld = ld+l1
                    i = iis
                    mesh = real(ld, kind=wp) * d_theta
                    product_1 = ZERO
                    do ii=3, ido, 2
                        i = i+2
                        product_1 = product_1 + ONE
                        theta = product_1 * mesh
                        trig_lookup_table(i-1) = cos(theta)
                        trig_lookup_table(i) = sin(theta)
                    end do
                    iis = iis+ido
                end do
                l1 = l2
            end do
        end if

    end subroutine precompute_factorization_and_trig_lookup_table

    subroutine forward_lower_utility_routine(m, n, c, mdimc, wa, fac)

        ! Dummy arguments
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: n
        real(wp),    intent(inout)  :: c(:,:)
        integer(ip), intent(in)     :: mdimc
        real(wp),    intent(in)     :: wa(:)
        real(wp),    intent(in)     :: fac(:)

        ! Local variables
        integer(ip) :: k1, n1, n2
        integer(ip) :: na, kh, nf, iip
        integer(ip) :: iw1, iw2, iw3, iw4, ido, idl1

        nf = int(fac(2), kind=ip)
        na = 1
        n2 = n
        iw1 = n

        block
            real(wp) :: ch(m,n)

            do k1=1, nf
                kh = nf-k1
                iip = int(fac(kh+3), kind=ip)
                n1 = n2/iip
                ido = n/n2
                idl1 = ido*n1
                iw1 = iw1-(iip-1)*ido
                na = 1-na

                select case (iip)
                    case (2)
                        if (na == 0) then
                            call forward_pass_2(m, ido, n1, c, mdimc, ch, m, wa(iw1:))
                        else
                            call forward_pass_2(m, ido, n1, ch, m, c, mdimc, wa(iw1:))
                        end if
                    case (3)
                        iw2 = iw1+ido
                        if (na == 0) then
                            call forward_pass_3(m, ido, n1, c, mdimc, ch, m, wa(iw1:), wa(iw2:))
                        else
                            call forward_pass_3(m, ido, n1, ch, m, c, mdimc, wa(iw1:), wa(iw2:))
                        end if
                    case(4)
                        iw2 = iw1+ido
                        iw3 = iw2+ido
                        if (na == 0) then
                            call forward_pass_4(m, ido, n1, c, mdimc, ch, m, wa(iw1:), wa(iw2:), wa(iw3:))
                        else
                            call forward_pass_4(m, ido, n1, ch, m, c, mdimc, wa(iw1:), wa(iw2:), wa(iw3:))
                        end if
                    case (5)
                        iw2 = iw1+ido
                        iw3 = iw2+ido
                        iw4 = iw3+ido
                        if (na == 0) then
                            call forward_pass_5(m, ido, n1, c, mdimc, ch, m, wa(iw1:), wa(iw2:), wa(iw3:), wa(iw4:))
                        else
                            call forward_pass_5(m, ido, n1, ch, m, c, mdimc, wa(iw1:), wa(iw2:), wa(iw3:), wa(iw4:))
                        end if
                    case default
                        if (ido == 1) na = 1-na
                        if (na == 0) then
                            call forward_pass_n(m, ido, iip, n1, idl1, c, c, c, mdimc, ch, ch, m, wa(iw1:))
                            na = 1
                        else
                            call forward_pass_n(m, ido, iip, n1, idl1, ch, ch, ch, m, c, c, mdimc, wa(iw1:))
                            na = 0
                        end if
                end select
                n2 = n1
            end do

            if (na /= 1) c(:m, :n) = ch
        end block

    end subroutine forward_lower_utility_routine

    subroutine forward_pass_2(mp, ido, l1, cc, mdimcc, ch, mdimch, wa1)

        ! Dummy arguments
        integer(ip), intent(in)     :: mp
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: l1
        real(wp),    intent(inout)  :: ch(mdimch, ido, 2, l1)
        integer(ip), intent(in)     :: mdimch
        real(wp),    intent(inout)  :: cc(mdimcc, ido, l1, 2)
        integer(ip), intent(in)     :: mdimcc
        real(wp),    intent(in)     :: wa1(:)

        ! Local variables
        integer(ip) :: i, k, m, ic, idp2

        ch(:mp, 1, 1, :) = cc(:mp, 1, :, 1)+cc(:mp, 1, :, 2)
        ch(:mp, ido, 2, :) = cc(:mp, 1, :, 1)-cc(:mp, 1, :, 2)

        if (ido < 2) then
            return
        else if (ido /= 2) then
            idp2 = ido+2
            do k=1, l1
                do i=3, ido, 2
                    ic = idp2-i
                    do m=1, mp
                        ch(m, i, 1, k) = &
                            cc(m, i, k, 1)+(wa1(i-2)*cc(m, i, k, 2)- &
                            wa1(i-1)*cc(m, i-1, k, 2))

                        ch(m, ic, 2, k) = &
                            (wa1(i-2)*cc(m, i, k, 2)-wa1(i-1)* &
                            cc(m, i-1, k, 2))-cc(m, i, k, 1)

                        ch(m, i-1, 1, k) = &
                            cc(m, i-1, k, 1)+(wa1(i-2)*cc(m, i-1, k, 2)+ &          !
                            wa1(i-1)*cc(m, i, k, 2))

                        ch(m, ic-1, 2, k) = &
                            cc(m, i-1, k, 1)-(wa1(i-2)*cc(m, i-1, k, 2)+ &         !
                            wa1(i-1)*cc(m, i, k, 2))
                    end do
                end do
            end do
            if (odd(ido)) return
        end if

        ch(:mp, 1, 2, :) = -cc(:mp, ido, :, 2)
        ch(:mp, ido, 1, :) = cc(:mp, ido, :, 1)

    end subroutine forward_pass_2

    subroutine forward_pass_3(mp, ido, l1, cc, mdimcc, ch, mdimch, wa1, wa2)

        ! Dummy arguments
        integer(ip), intent(in)     :: mp
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: l1
        real(wp),    intent(inout)  :: ch(mdimch, ido, 3, l1)
        integer(ip), intent(in)     :: mdimch
        real(wp),    intent(inout)  :: cc(mdimcc, ido, l1, 3)
        integer(ip), intent(in)     :: mdimcc
        real(wp),    intent(in)     :: wa1(:)
        real(wp),    intent(in)     :: wa2(:)

        ! Local variables
        integer(ip)         :: i, k, m, ic, idp2
        real(wp), parameter :: ARG=TWO_PI/3
        real(wp), parameter :: TAUR = cos(ARG)
        real(wp), parameter :: TAUI = sin(ARG)

        do k=1, l1
            do m=1, mp
                ch(m, 1, 1, k) = &
                    cc(m, 1, k, 1)+(cc(m, 1, k, 2)+cc(m, 1, k, 3))

                ch(m, 1, 3, k) = &
                    TAUI*(cc(m, 1, k, 3)-cc(m, 1, k, 2))

                ch(m, ido, 2, k) = &
                    cc(m, 1, k, 1)+TAUR* &
                    (cc(m, 1, k, 2)+cc(m, 1, k, 3))
            end do
        end do

        if (ido == 1) return

        idp2 = ido+2
        do k=1, l1
            do i=3, ido, 2
                ic = idp2-i
                do m=1, mp

                    ch(m, i-1, 1, k) = &
                        cc(m, i-1, k, 1)+((wa1(i-2)*cc(m, i-1, k, 2)+ &         !
                        wa1(i-1)*cc(m, i, k, 2))+(wa2(i-2)*cc(m, i-1, k, 3)+wa2(i-1)* &        !
                        cc(m, i, k, 3)))

                    ch(m, i, 1, k) = &
                        cc(m, i, k, 1)+((wa1(i-2)*cc(m, i, k, 2)-wa1(i-1)* &      !
                        cc(m, i-1, k, 2))+(wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)* &
                        cc(m, i-1, k, 3)))

                    ch(m, i-1, 3, k) = &
                        (cc(m, i-1, k, 1)+TAUR*((wa1(i-2)* &
                        cc(m, i-1, k, 2)+wa1(i-1)*cc(m, i, k, 2))+(wa2(i-2)* &
                        cc(m, i-1, k, 3)+wa2(i-1)*cc(m, i, k, 3))))+(TAUI*((wa1(i-2)* &        !
                        cc(m, i, k, 2)-wa1(i-1)*cc(m, i-1, k, 2))-(wa2(i-2)* &
                        cc(m, i, k, 3)-wa2(i-1)*cc(m, i-1, k, 3))))

                    ch(m, ic-1, 2, k) = &
                        (cc(m, i-1, k, 1)+TAUR*((wa1(i-2)* &
                        cc(m, i-1, k, 2)+wa1(i-1)*cc(m, i, k, 2))+(wa2(i-2)* &
                        cc(m, i-1, k, 3)+wa2(i-1)*cc(m, i, k, 3))))-(TAUI*((wa1(i-2)* &
                        cc(m, i, k, 2)-wa1(i-1)*cc(m, i-1, k, 2))-(wa2(i-2)* &
                        cc(m, i, k, 3)-wa2(i-1)*cc(m, i-1, k, 3))))

                    ch(m, i, 3, k) = &
                        (cc(m, i, k, 1)+TAUR*((wa1(i-2)*cc(m, i, k, 2)- &
                        wa1(i-1)*cc(m, i-1, k, 2))+(wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)* &
                        cc(m, i-1, k, 3))))+(TAUI*((wa2(i-2)*cc(m, i-1, k, 3)+wa2(i-1)* &
                        cc(m, i, k, 3))-(wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)* &
                        cc(m, i, k, 2))))

                    ch(m, ic, 2, k) = &
                        (TAUI*((wa2(i-2)*cc(m, i-1, k, 3)+wa2(i-1)* &
                        cc(m, i, k, 3))-(wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)* &
                        cc(m, i, k, 2))))-(cc(m, i, k, 1)+TAUR*((wa1(i-2)*cc(m, i, k, 2)- &
                        wa1(i-1)*cc(m, i-1, k, 2))+(wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)* &
                        cc(m, i-1, k, 3))))
                end do
            end do
        end do

    end subroutine forward_pass_3

    subroutine forward_pass_4(mp, ido, l1, cc, mdimcc, ch, mdimch, wa1, wa2, wa3)

        ! Dummy arguments
        integer(ip), intent(in)     :: mp
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: l1
        real(wp),    intent(inout)  :: cc(mdimcc, ido, l1, 4)
        integer(ip), intent(in)     :: mdimcc
        real(wp),    intent(inout)  :: ch(mdimch, ido, 4, l1)
        integer(ip), intent(in)     :: mdimch
        real(wp),    intent(in)     :: wa1(:)
        real(wp),    intent(in)     :: wa2(:)
        real(wp),    intent(in)     :: wa3(:)

        ! Local variables
        integer(ip)         :: i, k, m, ic, idp2
        real(wp), parameter :: HALF_SQRT_2 = SQRT_2/2

        do k=1, l1
            do m=1, mp
                ch(m, 1, 1, k) = &
                    (cc(m, 1, k, 2)+cc(m, 1, k, 4)) &
                    +(cc(m, 1, k, 1) + cc(m, 1, k, 3))
                ch(m, ido, 4, k) = &
                    (cc(m, 1, k, 1)+cc(m, 1, k, 3)) &
                    -(cc(m, 1, k, 2) + cc(m, 1, k, 4))
                ch(m, ido, 2, k) = &
                    cc(m, 1, k, 1)-cc(m, 1, k, 3)
                ch(m, 1, 3, k) = &
                    cc(m, 1, k, 4)-cc(m, 1, k, 2)
            end do
        end do

        if (ido < 2) then
            return
        else if (ido /= 2) then
            idp2 = ido+2
            do k=1, l1
                do i=3, ido, 2
                    ic = idp2-i
                    do m=1, mp

                        ch(m, i-1, 1, k) = &
                            ((wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)* &              !
                            cc(m, i, k, 2))+(wa3(i-2)*cc(m, i-1, k, 4)+wa3(i-1)* &
                            cc(m, i, k, 4)))+(cc(m, i-1, k, 1)+(wa2(i-2)*cc(m, i-1, k, 3)+ &          !
                            wa2(i-1)*cc(m, i, k, 3)))

                        ch(m, ic-1, 4, k) = &
                            (cc(m, i-1, k, 1)+(wa2(i-2)*cc(m, i-1, k, 3)+ &        !
                            wa2(i-1)*cc(m, i, k, 3)))-((wa1(i-2)*cc(m, i-1, k, 2)+ &               !
                            wa1(i-1)*cc(m, i, k, 2))+(wa3(i-2)*cc(m, i-1, k, 4)+ &
                            wa3(i-1)*cc(m, i, k, 4)))

                        ch(m, i, 1, k) = &
                            ((wa1(i-2)*cc(m, i, k, 2)-wa1(i-1)* &
                            cc(m, i-1, k, 2))+(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &
                            cc(m, i-1, k, 4)))+(cc(m, i, k, 1)+(wa2(i-2)*cc(m, i, k, 3)- &            !
                            wa2(i-1)*cc(m, i-1, k, 3)))

                        ch(m, ic, 4, k) = &
                            ((wa1(i-2)*cc(m, i, k, 2)-wa1(i-1)* &
                            cc(m, i-1, k, 2))+(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &
                            cc(m, i-1, k, 4)))-(cc(m, i, k, 1)+(wa2(i-2)*cc(m, i, k, 3)- &            !
                            wa2(i-1)*cc(m, i-1, k, 3)))

                        ch(m, i-1, 3, k) = &
                            ((wa1(i-2)*cc(m, i, k, 2)-wa1(i-1)* &
                            cc(m, i-1, k, 2))-(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &
                            cc(m, i-1, k, 4)))+(cc(m, i-1, k, 1)-(wa2(i-2)*cc(m, i-1, k, 3)+ &        !
                            wa2(i-1)*cc(m, i, k, 3)))

                        ch(m, ic-1, 2, k) = &
                            (cc(m, i-1, k, 1)-(wa2(i-2)*cc(m, i-1, k, 3)+ &        !
                            wa2(i-1)*cc(m, i, k, 3)))-((wa1(i-2)*cc(m, i, k, 2)-wa1(i-1)* &        !
                            cc(m, i-1, k, 2))-(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &
                            cc(m, i-1, k, 4)))

                        ch(m, i, 3, k) = &
                            ((wa3(i-2)*cc(m, i-1, k, 4)+wa3(i-1)* &
                            cc(m, i, k, 4))-(wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)* &
                            cc(m, i, k, 2)))+(cc(m, i, k, 1)-(wa2(i-2)*cc(m, i, k, 3)- &              !
                            wa2(i-1)*cc(m, i-1, k, 3)))

                        ch(m, ic, 2, k) = &
                            ((wa3(i-2)*cc(m, i-1, k, 4)+wa3(i-1)* &               !
                            cc(m, i, k, 4))-(wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)* &
                            cc(m, i, k, 2)))-(cc(m, i, k, 1)-(wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)* &     !
                            cc(m, i-1, k, 3)))
                    end do
                end do
            end do

            if (odd(ido)) return

        end if

        do k=1, l1
            do m=1, mp

                ch(m, ido, 1, k) = &
                    (HALF_SQRT_2*(cc(m, ido, k, 2)-cc(m, ido, k, 4)))+ &
                    cc(m, ido, k, 1)

                ch(m, ido, 3, k) = &
                    cc(m, ido, k, 1)-(HALF_SQRT_2*(cc(m, ido, k, 2)- &
                    cc(m, ido, k, 4)))

                ch(m, 1, 2, k) = &
                    (-HALF_SQRT_2*(cc(m, ido, k, 2)+cc(m, ido, k, 4)))- &
                    cc(m, ido, k, 3)

                ch(m, 1, 4, k) = &
                    (-HALF_SQRT_2*(cc(m, ido, k, 2)+cc(m, ido, k, 4)))+ &
                    cc(m, ido, k, 3)
            end do
        end do

    end subroutine forward_pass_4

    subroutine forward_pass_5(mp, ido, l1, cc, mdimcc, ch, mdimch, &
        wa1, wa2, wa3, wa4)

        ! Dummy arguments
        integer(ip), intent(in)     :: mp
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: l1
        real(wp),    intent(inout)  :: ch(mdimch, ido, 5, l1)
        integer(ip), intent(in)     :: mdimch
        real(wp),    intent(inout)  :: cc(mdimcc, ido, l1, 5)
        integer(ip), intent(in)     :: mdimcc
        real(wp),    intent(in)     :: wa1(:)
        real(wp),    intent(in)     :: wa2(:)
        real(wp),    intent(in)     :: wa3(:)
        real(wp),    intent(in)     :: wa4(:)

        ! Local variables
        integer(ip)         :: i, k, m, ic, idp2
        real(wp), parameter :: ARG = TWO_PI/5
        real(wp), parameter :: TR11=cos(ARG)
        real(wp), parameter :: TI11=sin(ARG)
        real(wp), parameter :: TR12=cos(TWO*ARG)
        real(wp), parameter :: TI12=sin(TWO*ARG)

        do k=1, l1
            do m=1, mp

                ch(m, 1, 1, k) = &
                    cc(m, 1, k, 1)+(cc(m, 1, k, 5)+cc(m, 1, k, 2))+ &
                    (cc(m, 1, k, 4)+cc(m, 1, k, 3))

                ch(m, ido, 2, k) = &
                    cc(m, 1, k, 1)+TR11*(cc(m, 1, k, 5)+cc(m, 1, k, 2))+ &
                    TR12*(cc(m, 1, k, 4)+cc(m, 1, k, 3))

                ch(m, 1, 3, k) = &
                    TI11*(cc(m, 1, k, 5)-cc(m, 1, k, 2))+TI12* &
                    (cc(m, 1, k, 4)-cc(m, 1, k, 3))

                ch(m, ido, 4, k) = &
                    cc(m, 1, k, 1)+TR12*(cc(m, 1, k, 5)+cc(m, 1, k, 2))+ &        !
                    TR11*(cc(m, 1, k, 4)+cc(m, 1, k, 3))

                ch(m, 1, 5, k) = &
                    TI12*(cc(m, 1, k, 5)-cc(m, 1, k, 2))-TI11* &
                    (cc(m, 1, k, 4)-cc(m, 1, k, 3))
            end do
        end do

        if (ido == 1) return

        idp2 = ido+2
        do k=1, l1
            do i=3, ido, 2
                ic = idp2-i
                do m=1, mp

                    ch(m, i-1, 1, k) = &
                        cc(m, i-1, k, 1)+((wa1(i-2)*cc(m, i-1, k, 2)+ &         !
                        wa1(i-1)*cc(m, i, k, 2))+(wa4(i-2)*cc(m, i-1, k, 5)+wa4(i-1)* &        !
                        cc(m, i, k, 5)))+((wa2(i-2)*cc(m, i-1, k, 3)+wa2(i-1)* &               !
                        cc(m, i, k, 3))+(wa3(i-2)*cc(m, i-1, k, 4)+wa3(i-1)*cc(m, i, k, 4)))      !

                    ch(m, i, 1, k) = &
                        cc(m, i, k, 1)+((wa1(i-2)*cc(m, i, k, 2)-wa1(i-1)* &      !
                        cc(m, i-1, k, 2))+(wa4(i-2)*cc(m, i, k, 5)-wa4(i-1)* &
                        cc(m, i-1, k, 5)))+((wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)* &               !
                        cc(m, i-1, k, 3))+(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &
                        cc(m, i-1, k, 4)))

                    ch(m, i-1, 3, k) = &
                        cc(m, i-1, k, 1)+TR11* &
                        ( wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)*cc(m, i, k, 2) &
                        +wa4(i-2)*cc(m, i-1, k, 5)+wa4(i-1)*cc(m, i, k, 5))+TR12* &            !
                        ( wa2(i-2)*cc(m, i-1, k, 3)+wa2(i-1)*cc(m, i, k, 3) &
                        +wa3(i-2)*cc(m, i-1, k, 4)+wa3(i-1)*cc(m, i, k, 4))+TI11* &            !
                        ( wa1(i-2)*cc(m, i, k, 2)-wa1(i-1)*cc(m, i-1, k, 2) &
                        -(wa4(i-2)*cc(m, i, k, 5)-wa4(i-1)*cc(m, i-1, k, 5)))+TI12* &          !
                        ( wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)*cc(m, i-1, k, 3) &
                        -(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)*cc(m, i-1, k, 4)))

                    ch(m, ic-1, 2, k) = &
                        cc(m, i-1, k, 1)+TR11* &
                        ( wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)*cc(m, i, k, 2) &
                        +wa4(i-2)*cc(m, i-1, k, 5)+wa4(i-1)*cc(m, i, k, 5))+TR12* &            !
                        ( wa2(i-2)*cc(m, i-1, k, 3)+wa2(i-1)*cc(m, i, k, 3) &
                        +wa3(i-2)*cc(m, i-1, k, 4)+wa3(i-1)*cc(m, i, k, 4))-(TI11* &            !
                        ( wa1(i-2)*cc(m, i, k, 2)-wa1(i-1)*cc(m, i-1, k, 2) &
                        -(wa4(i-2)*cc(m, i, k, 5)-wa4(i-1)*cc(m, i-1, k, 5)))+TI12* &          !
                        ( wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)*cc(m, i-1, k, 3) &
                        -(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)*cc(m, i-1, k, 4))))                 !

                    ch(m, i, 3, k) = &
                        (cc(m, i, k, 1)+TR11*((wa1(i-2)*cc(m, i, k, 2)- &         !
                        wa1(i-1)*cc(m, i-1, k, 2))+(wa4(i-2)*cc(m, i, k, 5)-wa4(i-1)* &        !
                        cc(m, i-1, k, 5)))+TR12*((wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)* &          !
                        cc(m, i-1, k, 3))+(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &
                        cc(m, i-1, k, 4))))+(TI11*((wa4(i-2)*cc(m, i-1, k, 5)+ &               !
                        wa4(i-1)*cc(m, i, k, 5))-(wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)* &        !
                        cc(m, i, k, 2)))+TI12*((wa3(i-2)*cc(m, i-1, k, 4)+wa3(i-1)* &          !
                        cc(m, i, k, 4))-(wa2(i-2)*cc(m, i-1, k, 3)+wa2(i-1)* &
                        cc(m, i, k, 3))))

                    ch(m, ic, 2, k) = &
                        (TI11*((wa4(i-2)*cc(m, i-1, k, 5)+wa4(i-1)* &         !
                        cc(m, i, k, 5))-(wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)* &
                        cc(m, i, k, 2)))+TI12*((wa3(i-2)*cc(m, i-1, k, 4)+wa3(i-1)* &          !
                        cc(m, i, k, 4))-(wa2(i-2)*cc(m, i-1, k, 3)+wa2(i-1)* &
                        cc(m, i, k, 3))))-(cc(m, i, k, 1)+TR11*((wa1(i-2)*cc(m, i, k, 2)- &       !
                        wa1(i-1)*cc(m, i-1, k, 2))+(wa4(i-2)*cc(m, i, k, 5)-wa4(i-1)* &        !
                        cc(m, i-1, k, 5)))+TR12*((wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)* &          !
                        cc(m, i-1, k, 3))+(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &
                        cc(m, i-1, k, 4))))

                    ch(m, i-1, 5, k) = &
                        (cc(m, i-1, k, 1)+TR12*((wa1(i-2)* &
                        cc(m, i-1, k, 2)+wa1(i-1)*cc(m, i, k, 2))+(wa4(i-2)* &
                        cc(m, i-1, k, 5)+wa4(i-1)*cc(m, i, k, 5)))+TR11*((wa2(i-2)* &          !
                        cc(m, i-1, k, 3)+wa2(i-1)*cc(m, i, k, 3))+(wa3(i-2)* &
                        cc(m, i-1, k, 4)+wa3(i-1)*cc(m, i, k, 4))))+(TI12*((wa1(i-2)* &        !
                        cc(m, i, k, 2)-wa1(i-1)*cc(m, i-1, k, 2))-(wa4(i-2)*cc(m, i, k, 5)- &     !
                        wa4(i-1)*cc(m, i-1, k, 5)))-TI11*((wa2(i-2)*cc(m, i, k, 3)- &          !
                        wa2(i-1)*cc(m, i-1, k, 3))-(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &        !
                        cc(m, i-1, k, 4))))

                    ch(m, ic-1, 4, k) = &
                        (cc(m, i-1, k, 1)+TR12*((wa1(i-2)* &
                        cc(m, i-1, k, 2)+wa1(i-1)*cc(m, i, k, 2))+(wa4(i-2)* &
                        cc(m, i-1, k, 5)+wa4(i-1)*cc(m, i, k, 5)))+TR11*((wa2(i-2)* &
                        cc(m, i-1, k, 3)+wa2(i-1)*cc(m, i, k, 3))+(wa3(i-2)* &
                        cc(m, i-1, k, 4)+wa3(i-1)*cc(m, i, k, 4))))-(TI12*((wa1(i-2)* &
                        cc(m, i, k, 2)-wa1(i-1)*cc(m, i-1, k, 2))-(wa4(i-2)*cc(m, i, k, 5)- &
                        wa4(i-1)*cc(m, i-1, k, 5)))-TI11*((wa2(i-2)*cc(m, i, k, 3)- &
                        wa2(i-1)*cc(m, i-1, k, 3))-(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &
                        cc(m, i-1, k, 4))))

                    ch(m, i, 5, k) = &
                        (cc(m, i, k, 1)+TR12*((wa1(i-2)*cc(m, i, k, 2)- &
                        wa1(i-1)*cc(m, i-1, k, 2))+(wa4(i-2)*cc(m, i, k, 5)-wa4(i-1)* &
                        cc(m, i-1, k, 5)))+TR11*((wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)* &
                        cc(m, i-1, k, 3))+(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &
                        cc(m, i-1, k, 4))))+(TI12*((wa4(i-2)*cc(m, i-1, k, 5)+ &
                        wa4(i-1)*cc(m, i, k, 5))-(wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)* &
                        cc(m, i, k, 2)))-TI11*((wa3(i-2)*cc(m, i-1, k, 4)+wa3(i-1)* &
                        cc(m, i, k, 4))-(wa2(i-2)*cc(m, i-1, k, 3)+wa2(i-1)* &
                        cc(m, i, k, 3))))

                    ch(m, ic, 4, k) = &
                        (TI12*((wa4(i-2)*cc(m, i-1, k, 5)+wa4(i-1)* &
                        cc(m, i, k, 5))-(wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)* &
                        cc(m, i, k, 2)))-TI11*((wa3(i-2)*cc(m, i-1, k, 4)+wa3(i-1)* &
                        cc(m, i, k, 4))-(wa2(i-2)*cc(m, i-1, k, 3)+wa2(i-1)* &
                        cc(m, i, k, 3))))-(cc(m, i, k, 1)+TR12*((wa1(i-2)*cc(m, i, k, 2)- &
                        wa1(i-1)*cc(m, i-1, k, 2))+(wa4(i-2)*cc(m, i, k, 5)-wa4(i-1)* &
                        cc(m, i-1, k, 5)))+TR11*((wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)* &
                        cc(m, i-1, k, 3))+(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &
                        cc(m, i-1, k, 4))))
                end do
            end do
        end do

    end subroutine forward_pass_5

    subroutine forward_pass_n(mp, ido, iip, l1, idl1, cc, c1, c2, mdimcc, &
        ch, ch2, mdimch, wa)

        ! Dummy arguments
        integer(ip), intent(in)     :: mp
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: iip
        integer(ip), intent(in)     :: l1
        integer(ip), intent(in)     :: idl1
        real(wp),    intent(inout)  :: cc(mdimcc, ido, iip, l1)
        real(wp),    intent(inout)  :: c1(mdimcc, ido, l1, iip)
        real(wp),    intent(inout)  :: c2(mdimcc, idl1, iip)
        integer(ip), intent(in)     :: mdimcc
        real(wp),    intent(inout)  :: ch(mdimch, ido, l1, iip)
        real(wp),    intent(inout)  :: ch2(mdimch, idl1, iip)
        integer(ip), intent(in)     :: mdimch
        real(wp),    intent(in)     :: wa(:)

        ! Local variables
        integer(ip) :: i, j, k, l, j2, ic, jc, lc, ik, iis, idij
        real(wp)    :: dc2, ai1, ai2, ar1, ar2, ds2
        real(wp)    :: ar1h, ar2h

        associate (&
            ipph => (iip+1)/2, &
            ipp2 => iip+2, &
            idp2 => ido+2, &
            nbd => (ido-1)/2, &
            arg => TWO_PI/iip &
            )
            associate (&
                dcp => cos(arg), &
                dsp => sin(arg) &
                )

                if (ido /= 1) then
                    ch2(:mp, :idl1, 1) = c2(:mp, :idl1, 1)
                    ch(:mp, 1, :l1, 2:iip) = c1(:mp, 1, :l1, 2:iip)

                    if (nbd <= l1) then
                        iis = -ido
                        do j=2, iip
                            iis = iis+ido
                            idij = iis
                            do i=3, ido, 2
                                ch(:mp, i-1, :l1, j) = wa(idij+1)*c1(:mp, i-1, :l1, j)+wa(idij+2) &            !
                                    *c1(:mp, i, :l1, j)
                                ch(:mp, i, :l1, j) = wa(idij+1)*c1(:mp, i, :l1, j)-wa(idij+2) &
                                    *c1(:mp, i-1, :l1, j)
                            end do
                        end do
                    else
                        iis = -ido
                        do j=2, iip
                            iis = iis+ido
                            do k=1, l1
                                idij = iis
                                do i=3, ido, 2
                                    idij = idij+2
                                    ch(:mp, i-1, k, j) = wa(idij-1)*c1(:mp, i-1, k, j)+wa(idij) &            !
                                        *c1(:mp, i, k, j)
                                    ch(:mp, i, k, j) = wa(idij-1)*c1(:mp, i, k, j)-wa(idij) &
                                        *c1(:mp, i-1, k, j)
                                end do
                            end do
                        end do
                    end if
                    if (l1 <= nbd) then
                        do j=2, ipph
                            jc = ipp2-j
                            do k=1, l1
                                do i=3, ido, 2
                                    c1(:mp, i-1, k, j) = ch(:mp, i-1, k, j)+ch(:mp, i-1, k, jc)
                                    c1(:mp, i-1, k, jc) = ch(:mp, i, k, j)-ch(:mp, i, k, jc)
                                    c1(:mp, i, k, j) = ch(:mp, i, k, j)+ch(:mp, i, k, jc)
                                    c1(:mp, i, k, jc) = ch(:mp, i-1, k, jc)-ch(:mp, i-1, k, j)
                                end do
                            end do
                        end do
                    else
                        do j=2, ipph
                            jc = ipp2-j
                            do i=3, ido, 2
                                c1(:mp, i-1, :l1, j) = ch(:mp, i-1, :l1, j)+ch(:mp, i-1, :l1, jc)
                                c1(:mp, i-1, :l1, jc) = ch(:mp, i, :l1, j)-ch(:mp, i, :l1, jc)
                                c1(:mp, i, :l1, j) = ch(:mp, i, :l1, j)+ch(:mp, i, :l1, jc)
                                c1(:mp, i, :l1, jc) = ch(:mp, i-1, :l1, jc)-ch(:mp, i-1, :l1, j)
                            end do
                        end do
                    end if
                else
                    c2(:mp, :idl1, 1) = ch2(:mp, :idl1, 1)
                end if

                do j=2, ipph
                    jc = ipp2-j
                    c1(:mp, 1, :l1, j) = ch(:mp, 1, :l1, j)+ch(:mp, 1, :l1, jc)
                    c1(:mp, 1, :l1, jc) = ch(:mp, 1, :l1, jc)-ch(:mp, 1, :l1, j)
                end do

                ar1 = ONE
                ai1 = ZERO
                do l=2, ipph
                    lc = ipp2-l
                    ar1h = dcp*ar1-dsp*ai1
                    ai1 = dcp*ai1+dsp*ar1
                    ar1 = ar1h
                    ch2(:mp, :idl1, l) = c2(:mp, :idl1, 1)+ar1*c2(:mp, :idl1, 2)
                    ch2(:mp, :idl1, lc) = ai1*c2(:mp, :idl1, iip)
                    dc2 = ar1
                    ds2 = ai1
                    ar2 = ar1
                    ai2 = ai1
                    do j=3, ipph
                        jc = ipp2-j
                        ar2h = dc2*ar2-ds2*ai2
                        ai2 = dc2*ai2+ds2*ar2
                        ar2 = ar2h
                        ch2(:mp, :idl1, l) = ch2(:mp, :idl1, l)+ar2*c2(:mp, :idl1, j)
                        ch2(:mp, :idl1, lc) = ch2(:mp, :idl1, lc)+ai2*c2(:mp, :idl1, jc)
                    end do
                end do

                do j=2, ipph
                    do ik=1, idl1
                        ch2(:mp, ik, 1) = ch2(:mp, ik, 1)+c2(:mp, ik, j)
                    end do
                end do

                if (ido >= l1) then
                    cc(:mp, :ido, 1, :l1) = ch(:mp, :ido, :l1, 1)
                else
                    cc(:mp, :ido, 1, :l1) = ch(:mp, :ido, :l1, 1)
                end if

                do j=2, ipph
                    jc = ipp2-j
                    j2 = 2*j
                    cc(:mp, ido, j2-2, :l1) = ch(:mp, 1, :l1, j)
                    cc(:mp, 1, j2-1, :l1) = ch(:mp, 1, :l1, jc)
                end do

                if (ido == 1) return

                if (l1 <= nbd) then
                    do j=2, ipph
                        jc = ipp2-j
                        j2 = 2*j
                        do k=1, l1
                            do i=3, ido, 2
                                ic = idp2-i
                                cc(:mp, i-1, j2-1, k) = ch(:mp, i-1, k, j)+ch(:mp, i-1, k, jc)                !
                                cc(:mp, ic-1, j2-2, k) = ch(:mp, i-1, k, j)-ch(:mp, i-1, k, jc)               !
                                cc(:mp, i, j2-1, k) = ch(:mp, i, k, j)+ch(:mp, i, k, jc)
                                cc(:mp, ic, j2-2, k) = ch(:mp, i, k, jc)-ch(:mp, i, k, j)
                            end do
                        end do
                    end do
                else
                    do j=2, ipph
                        jc = ipp2-j
                        j2 = 2*j
                        do i=3, ido, 2
                            ic = idp2-i
                            cc(:mp, i-1, j2-1, :l1) = ch(:mp, i-1, :l1, j)+ch(:mp, i-1, :l1, jc)                !
                            cc(:mp, ic-1, j2-2, :l1) = ch(:mp, i-1, :l1, j)-ch(:mp, i-1, :l1, jc)               !
                            cc(:mp, i, j2-1, :l1) = ch(:mp, i, :l1, j)+ch(:mp, i, :l1, jc)
                            cc(:mp, ic, j2-2, :l1) = ch(:mp, i, :l1, jc)-ch(:mp, i, :l1, j)
                        end do
                    end do
                end if
            end associate
        end associate

    end subroutine forward_pass_n

    subroutine backward_lower_utility_routine(m, n, c, mdimc, wa, fac)

        ! Dummy arguments
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: n
        real(wp),    intent(inout)  :: c(:,:)
        integer(ip), intent(in)     :: mdimc
        real(wp),    intent(in)     :: wa(:)
        real(wp),    intent(in)     :: fac(:)

        ! Local variables
        integer(ip) :: k1, n1, n2, na
        integer(ip) :: nf, iip, iw1, iw2, iw3, iw4, ido, idl1

        nf = int(fac(2), kind=ip)
        na = 0
        n1 = 1
        iw1 = 1

        block
            real(wp) :: ch(m,n)

            do k1=1, nf
                iip = int(fac(k1+2), kind=ip)
                n2 = iip*n1
                ido = n/n2
                idl1 = ido*n1

                select case (iip)
                    case (2)
                        if (na == 0) then
                            call backward_pass_2(m, ido, n1, c, mdimc, ch, m, wa(iw1:))
                        else
                            call backward_pass_2(m, ido, n1, ch, m, c, mdimc, wa(iw1:))
                        end if
                    case (3)
                        iw2 = iw1+ido
                        if (na == 0) then
                            call backward_pass_3(m, ido, n1, c, mdimc, ch, m, wa(iw1:), wa(iw2:))
                        else
                            call backward_pass_3(m, ido, n1, ch, m, c, mdimc, wa(iw1:), wa(iw2:))
                        end if
                    case (4)
                        iw2 = iw1+ido
                        iw3 = iw2+ido
                        if (na == 0) then
                            call backward_pass_4(m, ido, n1, c, mdimc, ch, m, wa(iw1:), wa(iw2:), wa(iw3:))
                        else
                            call backward_pass_4(m, ido, n1, ch, m, c, mdimc, wa(iw1:), wa(iw2:), wa(iw3:))
                        end if
                    case (5)
                        iw2 = iw1+ido
                        iw3 = iw2+ido
                        iw4 = iw3+ido
                        if (na == 0) then
                            call backward_pass_5(m, ido, n1, c, mdimc, ch, m, wa(iw1:), wa(iw2:), wa(iw3:), wa(iw4:))
                        else
                            call backward_pass_5(m, ido, n1, ch, m, c, mdimc, wa(iw1:), wa(iw2:), wa(iw3:), wa(iw4:))
                        end if
                    case default
                        if (na == 0) then
                            call backward_pass_n(m, ido, iip, n1, idl1, c, c, c, mdimc, ch, ch, m, wa(iw1:))
                        else
                            call backward_pass_n(m, ido, iip, n1, idl1, ch, ch, ch, m, c, c, mdimc, wa(iw1:))
                        end if
                        if (ido /= 1) na = 1-na
                end select
                na = 1-na
                n1 = n2
                iw1 = iw1+(iip-1)*ido
            end do

            if (na /= 0) c(:m, :n) = ch
        end block

    end subroutine backward_lower_utility_routine

    subroutine backward_pass_2(mp, ido, l1, cc, mdimcc, ch, mdimch, wa1)

        ! Dummy arguments
        integer(ip), intent(in)     :: mp
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: l1
        real(wp),    intent(inout)  :: cc(mdimcc, ido, 2, l1)
        integer(ip), intent(in)     :: mdimcc
        real(wp),    intent(out)    :: ch(mdimch, ido, l1, 2)
        integer(ip), intent(in)     :: mdimch
        real(wp),    intent(in)     :: wa1(:)

        ! Local variables
        integer(ip) :: i, k, ic, idp2

        ch(:mp, 1, :, 1) = cc(:mp, 1, 1, :)+cc(:mp, ido, 2, :)
        ch(:mp, 1, :, 2) = cc(:mp, 1, 1, :)-cc(:mp, ido, 2, :)

        if (ido < 2) then
            return
        else if (ido /= 2) then
            idp2 = ido+2
            do k=1, l1
                do i=3, ido, 2
                    ic = idp2-i
                    ch(:mp, i-1, k, 1) = &
                        cc(:mp, i-1, 1, k)+cc(:mp, ic-1, 2, k)

                    ch(:mp, i, k, 1) = &
                        cc(:mp, i, 1, k)-cc(:mp, ic, 2, k)

                    ch(:mp, i-1, k, 2) = &
                        wa1(i-2)*(cc(:mp, i-1, 1, k)-cc(:mp, ic-1, 2, k)) &         !
                        -wa1(i-1)*(cc(:mp, i, 1, k)+cc(:mp, ic, 2, k))

                    ch(:mp, i, k, 2) = &
                        wa1(i-2)*(cc(:mp, i, 1, k)+cc(:mp, ic, 2, k))+wa1(i-1) &      !
                        *(cc(:mp, i-1, 1, k)-cc(:mp, ic-1, 2, k))
                end do
            end do
            if (odd(ido)) return
        end if

        ch(:mp, ido, :, 1) = cc(:mp, ido, 1, :)+cc(:mp, ido, 1, :)
        ch(:mp, ido, :, 2) = -(cc(:mp, 1, 2, :)+cc(:mp, 1, 2, :))

    end subroutine backward_pass_2

    subroutine backward_pass_3(mp, ido, l1, cc, mdimcc, ch, mdimch, wa1, wa2)

        ! Dummy arguments
        integer(ip), intent(in)     :: mp
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: l1
        real(wp),    intent(inout)  :: cc(mdimcc, ido, 3, l1)
        integer(ip), intent(in)     :: mdimcc
        real(wp),    intent(out)    :: ch(mdimch, ido, l1, 3)
        integer(ip), intent(in)     :: mdimch
        real(wp),    intent(in)     :: wa1(:)
        real(wp),    intent(in)     :: wa2(:)

        ! Local variables
        integer(ip)         :: i, k, ic, idp2
        real(wp), parameter :: ARG = TWO_PI/3
        real(wp), parameter :: TAUR = cos(ARG)
        real(wp), parameter :: TAUI = sin(ARG)

        ch(:mp, 1, :l1, 1) = &
            cc(:mp, 1, 1, :l1)+TWO * cc(:mp, ido, 2, :l1)

        ch(:mp, 1, :l1, 2) = &
            cc(:mp, 1, 1, :l1)+(TWO * TAUR)*cc(:mp, ido, 2, :l1) &
            -(TWO *TAUI)*cc(:mp, 1, 3, :l1)

        ch(:mp, 1, :l1, 3) = &
            cc(:mp, 1, 1, :l1)+(TWO * TAUR)*cc(:mp, ido, 2, :l1) &
            +TWO *TAUI*cc(:mp, 1, 3, :l1)

        if (ido == 1) return

        idp2 = ido+2
        do k=1, l1
            do i=3, ido, 2
                ic = idp2-i
                ch(:mp, i-1, k, 1) = &
                    cc(:mp, i-1, 1, k)+(cc(:mp, i-1, 3, k)+cc(:mp, ic-1, 2, k))      !

                ch(:mp, i, k, 1) = &
                    cc(:mp, i, 1, k)+(cc(:mp, i, 3, k)-cc(:mp, ic, 2, k))              !

                ch(:mp, i-1, k, 2) = &
                    wa1(i-2)* &
                    ((cc(:mp, i-1, 1, k)+TAUR*(cc(:mp, i-1, 3, k)+cc(:mp, ic-1, 2, k)))- &
                    (TAUI*(cc(:mp, i, 3, k)+cc(:mp, ic, 2, k)))) &
                    -wa1(i-1)* &
                    ((cc(:mp, i, 1, k)+TAUR*(cc(:mp, i, 3, k)-cc(:mp, ic, 2, k)))+ &
                    (TAUI*(cc(:mp, i-1, 3, k)-cc(:mp, ic-1, 2, k))))

                ch(:mp, i, k, 2) = &
                    wa1(i-2)* &
                    ((cc(:mp, i, 1, k)+TAUR*(cc(:mp, i, 3, k)-cc(:mp, ic, 2, k)))+ &
                    (TAUI*(cc(:mp, i-1, 3, k)-cc(:mp, ic-1, 2, k)))) &
                    +wa1(i-1)* &
                    ((cc(:mp, i-1, 1, k)+TAUR*(cc(:mp, i-1, 3, k)+cc(:mp, ic-1, 2, k)))- &
                    (TAUI*(cc(:mp, i, 3, k)+cc(:mp, ic, 2, k))))

                ch(:mp, i-1, k, 3) = &
                    wa2(i-2)* &
                    ((cc(:mp, i-1, 1, k)+TAUR*(cc(:mp, i-1, 3, k)+cc(:mp, ic-1, 2, k)))+ &
                    (TAUI*(cc(:mp, i, 3, k)+cc(:mp, ic, 2, k)))) &
                    -wa2(i-1)* &
                    ((cc(:mp, i, 1, k)+TAUR*(cc(:mp, i, 3, k)-cc(:mp, ic, 2, k)))- &
                    (TAUI*(cc(:mp, i-1, 3, k)-cc(:mp, ic-1, 2, k))))

                ch(:mp, i, k, 3) = &
                    wa2(i-2)* &
                    ((cc(:mp, i, 1, k)+TAUR*(cc(:mp, i, 3, k)-cc(:mp, ic, 2, k)))- &
                    (TAUI*(cc(:mp, i-1, 3, k)-cc(:mp, ic-1, 2, k)))) &
                    +wa2(i-1)* &
                    ((cc(:mp, i-1, 1, k)+TAUR*(cc(:mp, i-1, 3, k)+cc(:mp, ic-1, 2, k)))+ &
                    (TAUI*(cc(:mp, i, 3, k)+cc(:mp, ic, 2, k))))
            end do
        end do

    end subroutine backward_pass_3

    subroutine backward_pass_4(mp, ido, l1, cc, mdimcc, ch, mdimch, wa1, wa2, wa3)

        ! Dummy arguments
        integer(ip), intent(in)     :: mp
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: l1
        real(wp),    intent(inout)  :: cc(mdimcc, ido, 4, l1)
        integer(ip), intent(in)     :: mdimcc
        real(wp),    intent(out)    :: ch(mdimch, ido, l1, 4)
        integer(ip), intent(in)     :: mdimch
        real(wp),    intent(in)     :: wa1(:)
        real(wp),    intent(in)     :: wa2(:)
        real(wp),    intent(in)     :: wa3(:)

        ! Local variables
        integer(ip)         :: i, k, ic, idp2
        real(wp), parameter :: SQRT2 = sqrt(TWO)


        ch(:mp, 1, :l1, 3) = (cc(:mp, 1, 1, :l1)+cc(:mp, ido, 4, :l1)) &
            -(cc(:mp, ido, 2, :l1)+cc(:mp, ido, 2, :l1))
        ch(:mp, 1, :l1, 1) = (cc(:mp, 1, 1, :l1)+cc(:mp, ido, 4, :l1)) &
            +(cc(:mp, ido, 2, :l1)+cc(:mp, ido, 2, :l1))
        ch(:mp, 1, :l1, 4) = (cc(:mp, 1, 1, :l1)-cc(:mp, ido, 4, :l1)) &
            +(cc(:mp, 1, 3, :l1)+cc(:mp, 1, 3, :l1))
        ch(:mp, 1, :l1, 2) = (cc(:mp, 1, 1, :l1)-cc(:mp, ido, 4, :l1)) &
            -(cc(:mp, 1, 3, :l1)+cc(:mp, 1, 3, :l1))

        if (ido < 2) then
            return
        else if (ido /= 2) then
            idp2 = ido+2
            do k=1, l1
                do i=3, ido, 2
                    ic = idp2-i
                    ch(:mp, i-1, k, 1) = &
                        (cc(:mp, i-1, 1, k)+cc(:mp, ic-1, 4, k)) &
                        +(cc(:mp, i-1, 3, k)+cc(:mp, ic-1, 2, k))

                    ch(:mp, i, k, 1) = &
                        (cc(:mp, i, 1, k)-cc(:mp, ic, 4, k)) &
                        +(cc(:mp, i, 3, k)-cc(:mp, ic, 2, k))

                    ch(:mp, i-1, k, 2) = &
                        wa1(i-2)*((cc(:mp, i-1, 1, k)-cc(:mp, ic-1, 4, k)) &          !
                        -(cc(:mp, i, 3, k)+cc(:mp, ic, 2, k)))-wa1(i-1) &
                        *((cc(:mp, i, 1, k)+cc(:mp, ic, 4, k))+(cc(:mp, i-1, 3, k)-cc(:mp, ic-1, 2, k)))      !

                    ch(:mp, i, k, 2)= &
                        wa1(i-2)*((cc(:mp, i, 1, k)+cc(:mp, ic, 4, k)) &
                        +(cc(:mp, i-1, 3, k)-cc(:mp, ic-1, 2, k)))+wa1(i-1) &
                        *((cc(:mp, i-1, 1, k)-cc(:mp, ic-1, 4, k))-(cc(:mp, i, 3, k)+cc(:mp, ic, 2, k)))      !

                    ch(:mp, i-1, k, 3) = &
                        wa2(i-2)*((cc(:mp, i-1, 1, k)+cc(:mp, ic-1, 4, k)) &          !
                        -(cc(:mp, i-1, 3, k)+cc(:mp, ic-1, 2, k)))-wa2(i-1) &
                        *((cc(:mp, i, 1, k)-cc(:mp, ic, 4, k))-(cc(:mp, i, 3, k)-cc(:mp, ic, 2, k)))          !

                    ch(:mp, i, k, 3) = &
                        wa2(i-2)*((cc(:mp, i, 1, k)-cc(:mp, ic, 4, k)) &
                        -(cc(:mp, i, 3, k)-cc(:mp, ic, 2, k)))+wa2(i-1) &
                        *((cc(:mp, i-1, 1, k)+cc(:mp, ic-1, 4, k))-(cc(:mp, i-1, 3, k) &
                        +cc(:mp, ic-1, 2, k)))

                    ch(:mp, i-1, k, 4) = &
                        wa3(i-2)*((cc(:mp, i-1, 1, k)-cc(:mp, ic-1, 4, k)) &          !
                        +(cc(:mp, i, 3, k)+cc(:mp, ic, 2, k)))-wa3(i-1) &
                        *((cc(:mp, i, 1, k)+cc(:mp, ic, 4, k))-(cc(:mp, i-1, 3, k)-cc(:mp, ic-1, 2, k)))       !

                    ch(:mp, i, k, 4) = &
                        wa3(i-2)*((cc(:mp, i, 1, k)+cc(:mp, ic, 4, k)) &
                        -(cc(:mp, i-1, 3, k)-cc(:mp, ic-1, 2, k)))+wa3(i-1) &
                        *((cc(:mp, i-1, 1, k)-cc(:mp, ic-1, 4, k))+(cc(:mp, i, 3, k)+cc(:mp, ic, 2, k)))      !
                end do
            end do
            if (odd(ido)) return
        end if

        ch(:mp, ido, :l1, 1) = &
            (cc(:mp, ido, 1, :l1)+cc(:mp, ido, 3, :l1)) &
            +(cc(:mp, ido, 1, :l1)+cc(:mp, ido, 3, :l1))

        ch(:mp, ido, :l1, 2) = &
            SQRT2*((cc(:mp, ido, 1, :l1)-cc(:mp, ido, 3, :l1)) &
            -(cc(:mp, 1, 2, :l1)+cc(:mp, 1, 4, :l1)))

        ch(:mp, ido, :l1, 3) = &
            (cc(:mp, 1, 4, :l1)-cc(:mp, 1, 2, :l1)) &
            +(cc(:mp, 1, 4, :l1)-cc(:mp, 1, 2, :l1))

        ch(:mp, ido, :l1, 4) = &
            -SQRT2*((cc(:mp, ido, 1, :l1)-cc(:mp, ido, 3, :l1)) &
            +(cc(:mp, 1, 2, :l1)+cc(:mp, 1, 4, :l1)))

    end subroutine backward_pass_4

    subroutine backward_pass_5(mp, ido, l1, cc, mdimcc, ch, mdimch, &
        wa1, wa2, wa3, wa4)

        ! Dummy arguments
        integer(ip), intent(in)     :: mp
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: l1
        real(wp),    intent(out)    :: ch(mdimch, ido, l1, 5)
        integer(ip), intent(in)     :: mdimch
        real(wp),    intent(inout)  :: cc(mdimcc, ido, 5, l1)
        integer(ip), intent(in)     :: mdimcc
        real(wp),    intent(in)     :: wa1(:)
        real(wp),    intent(in)     :: wa2(:)
        real(wp),    intent(in)     :: wa3(:)
        real(wp),    intent(in)     :: wa4(:)

        ! Local variables
        integer(ip)         :: i, k, ic, idp2
        real(wp), parameter :: ARG = TWO_PI/5
        real(wp), parameter :: TR11 = cos(ARG)
        real(wp), parameter :: TI11 = sin(ARG)
        real(wp), parameter :: TR12 = cos(TWO * ARG)
        real(wp), parameter :: TI12 = sin(TWO * ARG)

        do k=1, l1
            ch(:mp, 1, k, 1) = &
                cc(:mp, 1, 1, k)+TWO *cc(:mp, ido, 2, k)+TWO *cc(:mp, ido, 4, k)          !

            ch(:mp, 1, k, 2) = &
                (cc(:mp, 1, 1, k)+TR11*TWO *cc(:mp, ido, 2, k) &
                +TR12*TWO *cc(:mp, ido, 4, k))-(TI11*TWO *cc(:mp, 1, 3, k) &
                +TI12*TWO *cc(:mp, 1, 5, k))

            ch(:mp, 1, k, 3) = &
                (cc(:mp, 1, 1, k)+TR12*TWO *cc(:mp, ido, 2, k) &
                +TR11*TWO *cc(:mp, ido, 4, k))-(TI12*TWO *cc(:mp, 1, 3, k) &
                -TI11*TWO *cc(:mp, 1, 5, k))

            ch(:mp, 1, k, 4) = &
                (cc(:mp, 1, 1, k)+TR12*TWO *cc(:mp, ido, 2, k) &
                +TR11*TWO *cc(:mp, ido, 4, k))+(TI12*TWO *cc(:mp, 1, 3, k) &
                -TI11*TWO *cc(:mp, 1, 5, k))

            ch(:mp, 1, k, 5) = &
                (cc(:mp, 1, 1, k)+TR11*TWO *cc(:mp, ido, 2, k) &
                +TR12*TWO *cc(:mp, ido, 4, k))+(TI11*TWO *cc(:mp, 1, 3, k) &
                +TI12*TWO *cc(:mp, 1, 5, k))
        end do

        if (ido == 1) return

        idp2 = ido+2

        do k=1, l1
            do i=3, ido, 2
                ic = idp2-i
                ch(:mp, i-1, k, 1) = &
                    cc(:mp, i-1, 1, k)+(cc(:mp, i-1, 3, k)+cc(:mp, ic-1, 2, k)) &
                    +(cc(:mp, i-1, 5, k)+cc(:mp, ic-1, 4, k))

                ch(:mp, i, k, 1) = &
                    cc(:mp, i, 1, k)+(cc(:mp, i, 3, k)-cc(:mp, ic, 2, k)) &
                    +(cc(:mp, i, 5, k)-cc(:mp, ic, 4, k))

                ch(:mp, i-1, k, 2) = &
                    wa1(i-2)*((cc(:mp, i-1, 1, k)+TR11* &
                    (cc(:mp, i-1, 3, k)+cc(:mp, ic-1, 2, k))+TR12 &
                    *(cc(:mp, i-1, 5, k)+cc(:mp, ic-1, 4, k)))-(TI11*(cc(:mp, i, 3, k) &
                    +cc(:mp, ic, 2, k))+TI12*(cc(:mp, i, 5, k)+cc(:mp, ic, 4, k)))) &
                    -wa1(i-1)*((cc(:mp, i, 1, k)+TR11*(cc(:mp, i, 3, k)-cc(:mp, ic, 2, k)) &
                    +TR12*(cc(:mp, i, 5, k)-cc(:mp, ic, 4, k)))+(TI11*(cc(:mp, i-1, 3, k) &
                    -cc(:mp, ic-1, 2, k))+TI12*(cc(:mp, i-1, 5, k)-cc(:mp, ic-1, 4, k))))

                ch(:mp, i, k, 2) = &
                    wa1(i-2)*((cc(:mp, i, 1, k)+TR11*(cc(:mp, i, 3, k) &          !
                    -cc(:mp, ic, 2, k))+TR12*(cc(:mp, i, 5, k)-cc(:mp, ic, 4, k))) &
                    +(TI11*(cc(:mp, i-1, 3, k)-cc(:mp, ic-1, 2, k))+TI12 &
                    *(cc(:mp, i-1, 5, k)-cc(:mp, ic-1, 4, k))))+wa1(i-1) &
                    *((cc(:mp, i-1, 1, k)+TR11*(cc(:mp, i-1, 3, k) &
                    +cc(:mp, ic-1, 2, k))+TR12*(cc(:mp, i-1, 5, k)+cc(:mp, ic-1, 4, k))) &
                    -(TI11*(cc(:mp, i, 3, k)+cc(:mp, ic, 2, k))+TI12 &
                    *(cc(:mp, i, 5, k)+cc(:mp, ic, 4, k))))

                ch(:mp, i-1, k, 3) = &
                    wa2(i-2) &
                    *((cc(:mp, i-1, 1, k)+TR12*(cc(:mp, i-1, 3, k)+cc(:mp, ic-1, 2, k)) &
                    +TR11*(cc(:mp, i-1, 5, k)+cc(:mp, ic-1, 4, k)))-(TI12*(cc(:mp, i, 3, k) &
                    +cc(:mp, ic, 2, k))-TI11*(cc(:mp, i, 5, k)+cc(:mp, ic, 4, k)))) &
                    -wa2(i-1) &
                    *((cc(:mp, i, 1, k)+TR12*(cc(:mp, i, 3, k)- &
                    cc(:mp, ic, 2, k))+TR11*(cc(:mp, i, 5, k)-cc(:mp, ic, 4, k))) &
                    +(TI12*(cc(:mp, i-1, 3, k)-cc(:mp, ic-1, 2, k))-TI11 &
                    *(cc(:mp, i-1, 5, k)-cc(:mp, ic-1, 4, k))))

                ch(:mp, i, k, 3) = &
                    wa2(i-2) &
                    *((cc(:mp, i, 1, k)+TR12*(cc(:mp, i, 3, k)- &
                    cc(:mp, ic, 2, k))+TR11*(cc(:mp, i, 5, k)-cc(:mp, ic, 4, k))) &
                    +(TI12*(cc(:mp, i-1, 3, k)-cc(:mp, ic-1, 2, k))-TI11 &
                    *(cc(:mp, i-1, 5, k)-cc(:mp, ic-1, 4, k)))) &
                    +wa2(i-1) &
                    *((cc(:mp, i-1, 1, k)+TR12*(cc(:mp, i-1, 3, k)+cc(:mp, ic-1, 2, k)) &
                    +TR11*(cc(:mp, i-1, 5, k)+cc(:mp, ic-1, 4, k)))-(TI12*(cc(:mp, i, 3, k) &
                    +cc(:mp, ic, 2, k))-TI11*(cc(:mp, i, 5, k)+cc(:mp, ic, 4, k))))

                ch(:mp, i-1, k, 4) = &
                    wa3(i-2) &
                    *((cc(:mp, i-1, 1, k)+TR12*(cc(:mp, i-1, 3, k)+cc(:mp, ic-1, 2, k)) &
                    +TR11*(cc(:mp, i-1, 5, k)+cc(:mp, ic-1, 4, k)))+(TI12*(cc(:mp, i, 3, k) &
                    +cc(:mp, ic, 2, k))-TI11*(cc(:mp, i, 5, k)+cc(:mp, ic, 4, k)))) &
                    -wa3(i-1) &
                    *((cc(:mp, i, 1, k)+TR12*(cc(:mp, i, 3, k)- &
                    cc(:mp, ic, 2, k))+TR11*(cc(:mp, i, 5, k)-cc(:mp, ic, 4, k))) &
                    -(TI12*(cc(:mp, i-1, 3, k)-cc(:mp, ic-1, 2, k))-TI11 &
                    *(cc(:mp, i-1, 5, k)-cc(:mp, ic-1, 4, k))))

                ch(:mp, i, k, 4) = &
                    wa3(i-2) &
                    *((cc(:mp, i, 1, k)+TR12*(cc(:mp, i, 3, k)- &
                    cc(:mp, ic, 2, k))+TR11*(cc(:mp, i, 5, k)-cc(:mp, ic, 4, k))) &
                    -(TI12*(cc(:mp, i-1, 3, k)-cc(:mp, ic-1, 2, k))-TI11 &
                    *(cc(:mp, i-1, 5, k)-cc(:mp, ic-1, 4, k)))) &
                    +wa3(i-1) &
                    *((cc(:mp, i-1, 1, k)+TR12*(cc(:mp, i-1, 3, k)+cc(:mp, ic-1, 2, k)) &
                    +TR11*(cc(:mp, i-1, 5, k)+cc(:mp, ic-1, 4, k)))+(TI12*(cc(:mp, i, 3, k) &
                    +cc(:mp, ic, 2, k))-TI11*(cc(:mp, i, 5, k)+cc(:mp, ic, 4, k))))

                ch(:mp, i-1, k, 5) = &
                    wa4(i-2) &
                    *((cc(:mp, i-1, 1, k)+TR11*(cc(:mp, i-1, 3, k)+cc(:mp, ic-1, 2, k)) &
                    +TR12*(cc(:mp, i-1, 5, k)+cc(:mp, ic-1, 4, k)))+(TI11*(cc(:mp, i, 3, k) &
                    +cc(:mp, ic, 2, k))+TI12*(cc(:mp, i, 5, k)+cc(:mp, ic, 4, k)))) &
                    -wa4(i-1) &
                    *((cc(:mp, i, 1, k)+TR11*(cc(:mp, i, 3, k)-cc(:mp, ic, 2, k)) &
                    +TR12*(cc(:mp, i, 5, k)-cc(:mp, ic, 4, k)))-(TI11*(cc(:mp, i-1, 3, k) &
                    -cc(:mp, ic-1, 2, k))+TI12*(cc(:mp, i-1, 5, k)-cc(:mp, ic-1, 4, k))))

                ch(:mp, i, k, 5) = &
                    wa4(i-2) &
                    *((cc(:mp, i, 1, k)+TR11*(cc(:mp, i, 3, k)-cc(:mp, ic, 2, k)) &
                    +TR12*(cc(:mp, i, 5, k)-cc(:mp, ic, 4, k)))-(TI11*(cc(:mp, i-1, 3, k) &
                    -cc(:mp, ic-1, 2, k))+TI12*(cc(:mp, i-1, 5, k)-cc(:mp, ic-1, 4, k)))) &
                    +wa4(i-1) &
                    *((cc(:mp, i-1, 1, k)+TR11*(cc(:mp, i-1, 3, k)+cc(:mp, ic-1, 2, k)) &
                    +TR12*(cc(:mp, i-1, 5, k)+cc(:mp, ic-1, 4, k)))+(TI11*(cc(:mp, i, 3, k) &
                    +cc(:mp, ic, 2, k))+TI12*(cc(:mp, i, 5, k)+cc(:mp, ic, 4, k))))
            end do
        end do

    end subroutine backward_pass_5

    subroutine backward_pass_n(mp, ido, iip, l1, idl1, cc, c1, c2, mdimcc, &
        ch, ch2, mdimch, wa)

        ! Dummy arguments
        integer(ip), intent(in)     :: mp
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: iip
        integer(ip), intent(in)     :: l1
        integer(ip), intent(in)     :: idl1
        real(wp),    intent(inout)  :: cc(mdimcc, ido, iip, l1)
        real(wp),    intent(inout)  :: c1(mdimcc, ido, l1, iip)
        real(wp),    intent(inout)  :: c2(mdimcc, idl1, iip)
        integer(ip), intent(in)     :: mdimcc
        real(wp),    intent(inout)  :: ch(mdimch, ido, l1, iip)
        real(wp),    intent(inout)  :: ch2(mdimch, idl1, iip)
        integer(ip), intent(in)     :: mdimch
        real(wp),    intent(in)     :: wa(:)

        ! Local variables
        integer(ip)  :: i, j, k, l, j2, ic, jc, lc, iis, nbd
        integer(ip)  :: idp2, ipp2, idij, ipph
        real(wp)     :: dc2, ai1, ai2, ar1, ar2, ds2
        real(wp)     :: dcp, arg, dsp, ar1h, ar2h

        arg = TWO_PI/iip
        dcp = cos(arg)
        dsp = sin(arg)
        idp2 = ido+2
        nbd = (ido-1)/2
        ipp2 = iip+2
        ipph = (iip+1)/2

        ch(:mp, :, :, 1) = cc(:mp, :, 1, :)

        do j=2, ipph
            jc = ipp2-j
            j2 = 2*j
            ch(:mp, 1, :, j) = cc(:mp, ido, j2-2, :)+cc(:mp, ido, j2-2, :)
            ch(:mp, 1, :, jc) = cc(:mp, 1, j2-1, :)+cc(:mp, 1, j2-1, :)
        end do

        if (ido /= 1) then
            if (l1 <= nbd) then
                do j=2, ipph
                    jc = ipp2-j
                    do k=1, l1
                        do i=3, ido, 2
                            ic = idp2-i
                            ch(:mp, i-1, k, j) = cc(:mp, i-1, 2*j-1, k)+cc(:mp, ic-1, 2*j-2, k)
                            ch(:mp, i-1, k, jc) = cc(:mp, i-1, 2*j-1, k)-cc(:mp, ic-1, 2*j-2, k)
                            ch(:mp, i, k, j) = cc(:mp, i, 2*j-1, k)-cc(:mp, ic, 2*j-2, k)
                            ch(:mp, i, k, jc) = cc(:mp, i, 2*j-1, k)+cc(:mp, ic, 2*j-2, k)
                        end do
                    end do
                end do
            else
                do j=2, ipph
                    jc = ipp2-j
                    do i=3, ido, 2
                        ic = idp2-i
                        ch(:mp, i-1, :, j) = cc(:mp, i-1, 2*j-1, :)+cc(:mp, ic-1, 2*j-2, :)
                        ch(:mp, i-1, :, jc) = cc(:mp, i-1, 2*j-1, :)-cc(:mp, ic-1, 2*j-2, :)
                        ch(:mp, i, :, j) = cc(:mp, i, 2*j-1, :)-cc(:mp, ic, 2*j-2, :)
                        ch(:mp, i, :, jc) = cc(:mp, i, 2*j-1, :)+cc(:mp, ic, 2*j-2, :)
                    end do
                end do
            end if
        end if

        ar1 = ONE
        ai1 = ZERO
        do l=2, ipph
            lc = ipp2-l
            ar1h = dcp*ar1-dsp*ai1
            ai1 = dcp*ai1+dsp*ar1
            ar1 = ar1h
            c2(:mp, :idl1, l) = ch2(:mp, :idl1, 1)+ar1*ch2(:mp, :idl1, 2)
            c2(:mp, :idl1, lc) = ai1*ch2(:mp, :idl1, iip)
            dc2 = ar1
            ds2 = ai1
            ar2 = ar1
            ai2 = ai1
            do j=3, ipph
                jc = ipp2-j
                ar2h = dc2*ar2-ds2*ai2
                ai2 = dc2*ai2+ds2*ar2
                ar2 = ar2h
                c2(:mp, :idl1, l) = c2(:mp, :idl1, l)+ar2*ch2(:mp, :idl1, j)
                c2(:mp, :idl1, lc) = c2(:mp, :idl1, lc)+ai2*ch2(:mp, :idl1, jc)
            end do
        end do

        do j=2, ipph
            ch2(:mp, :idl1, 1) = ch2(:mp, :idl1, 1)+ch2(:mp, :idl1, j)
        end do

        do j=2, ipph
            jc = ipp2-j
            ch(:mp, 1, :l1, j) = c1(:mp, 1, :l1, j)-c1(:mp, 1, :l1, jc)
            ch(:mp, 1, :l1, jc) = c1(:mp, 1, :l1, j)+c1(:mp, 1, :l1, jc)
        end do

        if (ido /= 1) then
            if (l1 <= nbd) then
                do j=2, ipph
                    jc = ipp2-j
                    do k=1, l1
                        do i=3, ido, 2
                            ch(:mp, i-1, k, j) = c1(:mp, i-1, k, j)-c1(:mp, i, k, jc)
                            ch(:mp, i-1, k, jc) = c1(:mp, i-1, k, j)+c1(:mp, i, k, jc)
                            ch(:mp, i, k, j) = c1(:mp, i, k, j)+c1(:mp, i-1, k, jc)
                            ch(:mp, i, k, jc) = c1(:mp, i, k, j)-c1(:mp, i-1, k, jc)
                        end do
                    end do
                end do
            else
                do j=2, ipph
                    jc = ipp2-j
                    do i=3, ido, 2
                        ch(:mp, i-1, :l1, j) = c1(:mp, i-1, :l1, j)-c1(:mp, i, :l1, jc)
                        ch(:mp, i-1, :l1, jc) = c1(:mp, i-1, :l1, j)+c1(:mp, i, :l1, jc)
                        ch(:mp, i, :l1, j) = c1(:mp, i, :l1, j)+c1(:mp, i-1, :l1, jc)
                        ch(:mp, i, :l1, jc) = c1(:mp, i, :l1, j)-c1(:mp, i-1, :l1, jc)
                    end do
                end do
            end if
        end if

        if (ido /= 1) then
            c2(:mp, :idl1, 1) = ch2(:mp, :idl1, 1)
            c1(:mp, 1, :l1, 2:iip) = ch(:mp, 1, :l1, 2:iip)
            if (nbd <= l1) then
                iis = -ido
                do j=2, iip
                    iis = iis+ido
                    idij = iis
                    do i=3, ido, 2
                        idij = idij+2
                        do k=1, l1
                            c1(:mp, i-1, k, j) = &
                                wa(idij-1) * ch(:mp, i-1, k, j) &
                                - wa(idij) * ch(:mp, i, k, j)
                            c1(:mp, i, k, j) = &
                                wa(idij-1) * ch(:mp, i, k, j) &
                                + wa(idij) * ch(:mp, i-1, k, j)
                        end do
                    end do
                end do
            else
                iis = -ido
                do j=2, iip
                    iis = iis+ido
                    do k=1, l1
                        idij = iis
                        do i=3, ido, 2
                            idij = idij+2
                            c1(:mp, i-1, k, j) = &
                                wa(idij-1) * ch(:mp, i-1, k, j) &
                                - wa(idij) * ch(:mp, i, k, j)
                            c1(:mp, i, k, j) = &
                                wa(idij-1) * ch(:mp, i, k, j) &
                                + wa(idij) * ch(:mp, i-1, k, j)
                        end do
                    end do
                end do
            end if
        end if

    end subroutine backward_pass_n

end module type_RealPeriodicFastFourierTransform
