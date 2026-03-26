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
submodule(gradient_routines) invert_gradient_gaussian_grid

contains

    !     subroutine igradgc(nlat, nlon, isym, nt, sf, isf, jsf, br, bi, mdb, ndb, &
    !                        wshsgc, ierror)
    !
    !     let br, bi, cr, ci be the vector spherical harmonic coefficients
    !     precomputed by vhagc for a vector field (v, w).  let (v', w') be
    !     the irrotational component of (v, w) (i.e., (v', w') is generated
    !     by assuming cr, ci are zero and synthesizing br, bi with vhsgs).
    !     then subroutine igradgc computes a scalar field sf such that
    !
    !            gradient(sf) = (v', w').
    !
    !     i.e.,
    !
    !            v'(i, j) = d(sf(i, j))/dtheta          (colatitudinal component of
    !                                                 the gradient)
    !     and
    !
    !            w'(i, j) = 1/sint*d(sf(i, j))/dlambda  (east longitudinal component
    !                                                 of the gradient)
    !
    !     at the gaussian colatitude theta(i) (see nlat as input parameter)
    !     and longitude lambda(j) = (j-1)*2*pi/nlon where sint = sin(theta(i)).
    !
    !     note:  for an irrotational vector field (v, w), subroutine igradgc
    !     computes a scalar field whose gradient is (v, w).  in ay case,
    !     subroutine igradgc "inverts" the gradient subroutine gradgc.
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
    !            nlon = 72 for a five degree grid. nlon must be greater than
    !            3.  the efficiency of the computation is improved when nlon
    !            is a product of small prime numbers.
    !
    !
    !     isym   a parameter which determines whether the scalar field sf is
    !            computed on the full or half sphere as follows:
    !
    !      = 0
    !
    !            the symmetries/antsymmetries described in isym=1, 2 below
    !            do not exist in (v, w) about the equator.  in this case sf
    !            is neither symmetric nor antisymmetric about the equator.
    !            sf is computed on the entire sphere.  i.e., in the array
    !            sf(i, j) for i=1, ..., nlat and  j=1, ..., nlon
    !
    !      = 1
    !
    !            w is antisymmetric and v is symmetric about the equator.
    !            in this case sf is antisymmetyric about the equator and
    !            is computed for the northern hemisphere only.  i.e.,
    !            if nlat is odd sf is computed in the array sf(i, j) for
    !            i=1, ..., (nlat + 1)/2 and for j=1, ..., nlon.  if nlat is even
    !            sf is computed in the array sf(i, j) for i=1, ..., nlat/2
    !            and j=1, ..., nlon.
    !
    !      = 2
    !
    !            w is symmetric and v is antisymmetric about the equator.
    !            in this case sf is symmetyric about the equator and
    !            is computed for the northern hemisphere only.  i.e.,
    !            if nlat is odd sf is computed in the array sf(i, j) for
    !            i=1, ..., (nlat + 1)/2 and for j=1, ..., nlon.  if nlat is even
    !            sf is computed in the array sf(i, j) for i=1, ..., nlat/2
    !            and j=1, ..., nlon.
    !
    !
    !     nt     nt is the number of scalar and vector fields.  some
    !            computational efficiency is obtained for multiple fields.
    !            the arrays br, bi, and sf can be three dimensional corresponding
    !            to an indexed multiple vector field (v, w).  in this case,
    !            multiple scalar synthesis will be performed to compute each
    !            scalar field.  the third index for br, bi, and sf is the synthesis
    !            index which assumes the values k = 1, ..., nt.  for a single
    !            synthesis set nt = 1.  the description of the remaining
    !            parameters is simplified by assuming that nt=1 or that br, bi,
    !            and sf are two dimensional arrays.
    !
    !     isf    the first dimension of the array sf as it appears in
    !            the program that calls igradgc. if isym = 0 then isf
    !            must be at least nlat.  if isym = 1 or 2 and nlat is
    !            even then isf must be at least nlat/2. if isym = 1 or 2
    !            and nlat is odd then isf must be at least (nlat + 1)/2.
    !
    !     jsf    the second dimension of the array sf as it appears in
    !            the program that calls igradgc. jsf must be at least nlon.
    !
    !     br, bi  two or three dimensional arrays (see input parameter nt)
    !            that contain vector spherical harmonic coefficients
    !            of the vector field (v, w) as computed by subroutine vhagc.
    !     ***    br, bi must be computed by vhagc prior to calling igradgc.
    !
    !     mdb    the first dimension of the arrays br and bi as it appears in
    !            the program that calls igradgc (and vhagc). mdb must be at
    !            least min(nlat, nlon/2) if nlon is even or at least
    !            min(nlat, (nlon + 1)/2) if nlon is odd.
    !
    !     ndb    the second dimension of the arrays br and bi as it appears in
    !            the program that calls igradgc (and vhagc). ndb must be at
    !            least nlat.
    !
    !
    !  wshsgc    an array which must be initialized by subroutine shsgci.
    !            once initialized,
    !            wshsgc can be used repeatedly by igradgc as long as nlon
    !            and nlat remain unchanged.  wshsgc must not be altered
    !            between calls of igradgc.
    !
    !
    !  lshsgc    the dimension of the array wshsgc as it appears in the
    !            program that calls igradgc. define
    !
    !               l1 = min(nlat, (nlon+2)/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd.
    !
    !
    !            then lshsgc must be at least
    !
    !               nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15
    !
    !     output parameters
    !
    !
    !     sf    a two or three dimensional array (see input parameter nt) that
    !           contain a scalar field whose gradient is the irrotational
    !           component of the vector field (v, w).  the vector spherical
    !           harmonic coefficients br, bi were precomputed by subroutine
    !           vhagc.  sf(i, j) is given at the gaussian colatitude theta(i)
    !           and longitude lambda(j) = (j-1)*2*pi/nlon.  the index ranges
    !           are defined at input parameter isym.
    !
    !
    !  ierror   = 0  no errors
    !           = 1  error in the specification of nlat
    !           = 2  error in the specification of nlon
    !           = 3  error in the specification of isym
    !           = 4  error in the specification of nt
    !           = 5  error in the specification of isf
    !           = 6  error in the specification of jsf
    !           = 7  error in the specification of mdb
    !           = 8  error in the specification of ndb
    !           = 9  error in the specification of lshsgc
    !
    module subroutine igradgc(nlat, nlon, isym, nt, sf, isf, jsf, br, bi, mdb, ndb, &
        wshsgc, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: isym
        integer(ip), intent(in)  :: nt
        real(wp),    intent(out) :: sf(isf, jsf, nt)
        integer(ip), intent(in)  :: isf
        integer(ip), intent(in)  :: jsf
        real(wp),    intent(in)  :: br(mdb, ndb, nt)
        real(wp),    intent(in)  :: bi(mdb, ndb, nt)
        integer(ip), intent(in)  :: mdb
        integer(ip), intent(in)  :: ndb
        real(wp),    intent(in)  :: wshsgc(:)
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip) :: imid, n1, n2, mmax
        integer(ip) :: required_wavetable_size

        ! Check calling arguments
        ierror = 1
        if (nlat < 3) return
        ierror = 2
        if (nlon < 4) return
        ierror = 3
        if (isym < 0 .or. isym > 2) return
        ierror = 4
        if (nt < 0) return
        ierror = 5
        imid = (nlat + 1)/2
        if ((isym == 0 .and. isf<nlat) .or. &
            (isym /= 0 .and. isf<imid)) return
        ierror = 6
        if (jsf < nlon) return
        ierror = 7
        mmax = min(nlat, (nlon+2)/2)
        if (mdb < min(nlat, (nlon + 1)/2)) return
        ierror = 8
        if (ndb < nlat) return
        ierror = 9
        !
        !     verify saved workspace length
        !
        n2 = (nlat+mod(nlat, 2))/2
        n1 = min((nlon+2)/2, nlat)
        required_wavetable_size = nlat*(2*n2+3*n1-2)+3*n1*(1-n1)/2+nlon+15
        if (size(wshsgc) < required_wavetable_size) return
        ierror = 0

        call invert_gradient_lower_utility_routine(nlat, nlon, isym, nt, sf, &
            br, bi, wshsgc, shsgc, ierror)

    end subroutine igradgc

end submodule invert_gradient_gaussian_grid
