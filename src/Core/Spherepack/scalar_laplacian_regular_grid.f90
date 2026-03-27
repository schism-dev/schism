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
submodule(scalar_laplacian_routines) scalar_laplacian_regular_grid

contains
    !     subroutine slapec(nlat, nlon, isym, nt, slap, ids, jds, a, b, mdab, ndab,
    !                       wshsec, ierror)
    !
    !
    !     given the scalar spherical harmonic coefficients a and b, precomputed
    !     by subroutine shaec for a scalar field sf, subroutine slapec computes
    !     the laplacian of sf in the scalar array slap.  slap(i, j) is the
    !     laplacian of sf at the colatitude
    !
    !         theta(i) = (i-1)*pi/(nlat-1)
    !
    !     and east longitude
    !
    !         lambda(j) = (j-1)*2*pi/nlon
    !
    !     on the sphere.  i.e.
    !
    !         slap(i, j) =
    !
    !                  2                2
    !         [1/sint*d (sf(i, j)/dlambda + d(sint*d(sf(i, j))/dtheta)/dtheta]/sint
    !
    !
    !     where sint = sin(theta(i)).  the scalar laplacian in slap has the
    !     same symmetry or absence of symmetry about the equator as the scalar
    !     field sf.  the input parameters isym, nt, mdab, ndab must have the
    !     same values used by shaec to compute a and b for sf. the associated
    !     legendre functions are recomputed rather than stored as they are
    !     in subroutine slapes.

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
    !     nlon   the number of distinct longitude points.  nlon determines
    !            the grid increment in longitude as 2*pi/nlon. for example
    !            nlon = 72 for a five degree grid. nlon must be greater
    !            than zero. the axisymmetric case corresponds to nlon=1.
    !            the efficiency of the computation is improved when nlon
    !            is a product of small prime numbers.
    !
    !     isym   this parameter should have the same value input to subroutine
    !            shaec to compute the coefficients a and b for the scalar field
    !            sf.  isym is set as follows:
    !
    !            = 0  no symmetries exist in sf about the equator. scalar
    !                 synthesis is used to compute slap on the entire sphere.
    !                 i.e., in the array slap(i, j) for i=1, ..., nlat and
    !                 j=1, ..., nlon.
    !
    !           = 1  sf and slap are antisymmetric about the equator. the
    !                synthesis used to compute slap is performed on the
    !                northern hemisphere only.  if nlat is odd, slap(i, j) is
    !                computed for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon.  if
    !                nlat is even, slap(i, j) is computed for i=1, ..., nlat/2
    !                and j=1, ..., nlon.
    !
    !
    !           = 2  sf and slap are symmetric about the equator. the
    !                synthesis used to compute slap is performed on the
    !                northern hemisphere only.  if nlat is odd, slap(i, j) is
    !                computed for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon.  if
    !                nlat is even, slap(i, j) is computed for i=1, ..., nlat/2
    !                and j=1, ..., nlon.
    !
    !
    !     nt     the number of analyses.  in the program that calls slapec
    !            the arrays slap, a, and b can be three dimensional in which
    !            case multiple synthesis will be performed.  the third index
    !            is the synthesis index which assumes the values k=1, ..., nt.
    !            for a single analysis set nt=1. the description of the
    !            remaining parameters is simplified by assuming that nt=1
    !            or that all the arrays are two dimensional.
    !
    !   ids      the first dimension of the array slap as it appears in the
    !            program that calls slapec.  if isym = 0 then ids must be at
    !            least nlat.  if isym > 0 and nlat is even then ids must be
    !            at least nlat/2. if isym > 0 and nlat is odd then ids must
    !            be at least (nlat + 1)/2.
    !
    !   jds      the second dimension of the array slap as it appears in the
    !            program that calls slapec. jds must be at least nlon.
    !
    !
    !   a, b      two or three dimensional arrays (see input parameter nt)
    !            that contain scalar spherical harmonic coefficients
    !            of the scalar field sf as computed by subroutine shaec.
    !     ***    a, b must be computed by shaec prior to calling slapec.
    !
    !
    !    mdab    the first dimension of the arrays a and b as it appears
    !            in the program that calls slapec.  mdab must be at
    !            least min(nlat, (nlon+2)/2) if nlon is even or at least
    !            min(nlat, (nlon + 1)/2) if nlon is odd.
    !
    !    ndab    the second dimension of the arrays a and b as it appears
    !            in the program that calls slapec. ndbc must be at least
    !            least nlat.
    !
    !            mdab, ndab should have the same values input to shaec to
    !            compute the coefficients a and b.
    !
    !
    !    wshsec  an array which must be initialized by subroutine shseci
    !            before calling slapec.  once initialized, wshsec
    !            can be used repeatedly by slapec as long as nlat and nlon
    !            remain unchanged.  wshsec must not be altered between calls
    !            of slapec.
    !
    !    lshsec  the dimension of the array wshsec as it appears in the
    !            program that calls slapec.  let
    !
    !               l1 = min(nlat, (nlon+2)/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lshsec must be greater than or equal to
    !
    !               2*nlat*l2+3*((l1-2)*(2*nlat-l1-1))/2+nlon+15
    !
    !     output parameters
    !
    !
    !    slap    a two or three dimensional arrays (see input parameter nt) that
    !            contain the scalar laplacian of the scalar field sf.  slap(i, j)
    !            is the scalar laplacian at the colatitude
    !
    !                 theta(i) = (i-1)*pi/(nlat-1)
    !
    !            and longitude
    !
    !                 lambda(j) = (j-1)*2*pi/nlon
    !
    !            for i=1, ..., nlat and j=1, ..., nlon.
    !
    !
    !  ierror    a parameter which flags errors in input parameters as follows:
    !
    !            = 0  no errors detected
    !            = 1  error in the specification of nlat
    !            = 2  error in the specification of nlon
    !            = 3  error in the specification of isym
    !            = 4  error in the specification of nt
    !            = 5  error in the specification of ids
    !            = 6  error in the specification of jds
    !            = 7  error in the specification of mdab
    !            = 8  error in the specification of ndab
    !            = 9  error in the specification of lshsec
    !
    module subroutine slapec(nlat, nlon, isym, nt, slap, ids, jds, a, b, mdab, ndab, &
        wshsec, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: isym
        integer(ip), intent(in)  :: nt
        real(wp),    intent(out) :: slap(ids, jds, nt)
        integer(ip), intent(in)  :: ids
        integer(ip), intent(in)  :: jds
        real(wp),    intent(in)  :: a(mdab, ndab, nt)
        real(wp),    intent(in)  :: b(mdab, ndab, nt)
        integer(ip), intent(in)  :: mdab
        integer(ip), intent(in)  :: ndab
        real(wp),    intent(in)  :: wshsec(:)
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip)                  :: required_wavetable_size
        type(ScalarSynthesisUtility) :: util

        ! Check input arguments
        required_wavetable_size = util%get_lshsec(nlat, nlon)

        call util%check_scalar_transform_inputs(isym, ids, jds, &
            mdab, ndab, nlat, nlon, nt, required_wavetable_size, &
            wshsec, ierror)

        ! Check error flag
        if (ierror /= 0) return

        call scalar_laplacian_lower_utility_routine(nlat, nlon, isym, nt, slap, &
            a, b, wshsec, shsec, ierror)

    end subroutine slapec

end submodule scalar_laplacian_regular_grid
