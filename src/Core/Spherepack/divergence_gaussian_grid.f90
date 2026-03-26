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
submodule(divergence_routines) divergence_gaussian_grid

contains

    !     subroutine divgc(nlat, nlon, isym, nt, dv, idv, jdv, br, bi, mdb, ndb, &
    !                      wshsgc, ierror)
    !
    !     given the vector spherical harmonic coefficients br and bi, precomputed
    !     by subroutine vhagc for a vector field (v, w), subroutine divgc
    !     computes the divergence of the vector field in the scalar array dv.
    !     dv(i, j) is the divergence at the gaussian colatitude point theta(i)
    !     (see nlat as input parameter) and east longitude
    !
    !            lambda(j) = (j-1)*2*pi/nlon
    !
    !     on the sphere.  i.e.
    !
    !            dv(i, j) = 1/sint*[ d(sint*v(i, j))/dtheta + d(w(i, j))/dlambda ]
    !
    !     where sint = sin(theta(i)).  w is the east longitudinal and v
    !     is the colatitudinal component of the vector field from which
    !     br, bi were precomputed.  required associated legendre polynomials
    !     are recomputed rather than stored as they are in subroutine divgs.
    !
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
    !
    !     isym   a parameter which determines whether the divergence is
    !            computed on the full or half sphere as follows:
    !
    !      = 0
    !
    !            the symmetries/antsymmetries described in isym=1, 2 below
    !            do not exist in (v, w) about the equator.  in this case the
    !            divergence is neither symmetric nor antisymmetric about
    !            the equator.  the divergence is computed on the entire
    !            sphere.  i.e., in the array dv(i, j) for i=1, ..., nlat and
    !            j=1, ..., nlon.
    !
    !      = 1
    !
    !            w is antisymmetric and v is symmetric about the equator.
    !            in this case the divergence is antisymmetyric about
    !            the equator and is computed for the northern hemisphere
    !            only.  i.e., if nlat is odd the divergence is computed
    !            in the array dv(i, j) for i=1, ..., (nlat + 1)/2 and for
    !            j=1, ..., nlon.  if nlat is even the divergence is computed
    !            in the array dv(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !      = 2
    !            w is symmetric and v is antisymmetric about the equator
    !            in this case the divergence is symmetyric about the
    !            equator and is computed for the northern hemisphere
    !            only.  i.e., if nlat is odd the divergence is computed
    !            in the array dv(i, j) for i=1, ..., (nlat + 1)/2 and for
    !            j=1, ..., nlon.  if nlat is even the divergence is computed
    !            in the array dv(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !
    !     nt     nt is the number of scalar and vector fields.  some
    !            computational efficiency is obtained for multiple fields.
    !            in the program that calls divgc, the arrays br, bi, and dv
    !            can be three dimensional corresponding to an indexed multiple
    !            vector field.  in this case multiple scalar synthesis will
    !            be performed to compute the divergence for each field.  the
    !            third index is the synthesis index which assumes the values
    !            k=1, ..., nt.  for a single synthesis set nt = 1.  the
    !            description of the remaining parameters is simplified by
    !            assuming that nt=1 or that all the arrays are two dimensional.
    !
    !     idv    the first dimension of the array dv as it appears in
    !            the program that calls divgc. if isym = 0 then idv
    !            must be at least nlat.  if isym = 1 or 2 and nlat is
    !            even then idv must be at least nlat/2. if isym = 1 or 2
    !            and nlat is odd then idv must be at least (nlat + 1)/2.
    !
    !     jdv    the second dimension of the array dv as it appears in
    !            the program that calls divgc. jdv must be at least nlon.
    !
    !     br, bi  two or three dimensional arrays (see input parameter nt)
    !            that contain vector spherical harmonic coefficients
    !            of the vector field (v, w) as computed by subroutine vhagc.
    !     ***    br and bi must be computed by vhagc prior to calling
    !            divgc.
    !
    !     mdb    the first dimension of the arrays br and bi as it
    !            appears in the program that calls divgc. mdb must be at
    !            least min(nlat, nlon/2) if nlon is even or at least
    !            min(nlat, (nlon + 1)/2) if nlon is odd.
    !
    !     ndb    the second dimension of the arrays br and bi as it
    !            appears in the program that calls divgc. ndb must be at
    !            least nlat.
    !
    !
    !     wshsgc an array which must be initialized by subroutine shsgci
    !            once initialized, wshsgc can be used repeatedly by divgc
    !            as long as nlon and nlat remain unchanged.  wshsgc must
    !            not be altered between calls of divgc.
    !
    !
    !     lshsgc the dimension of the array wshsgc as it appears in the
    !            program that calls divgc. define
    !
    !               l1 = min(nlat, (nlon+2)/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lshsgc must be at least
    !
    !               nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15
    !
    !     output parameters
    !
    !
    !    dv     a two or three dimensional array (see input parameter nt)
    !           that contains the divergence of the vector field (v, w)
    !           whose coefficients br, bi where computed by subroutine
    !           vhagc.  dv(i, j) is the divergence at the gaussian colatitude
    !           point theta(i) and longitude point lambda(j) = (j-1)*2*pi/nlon.
    !           the index ranges are defined above at the input parameter
    !           isym.
    !
    !
    !    ierror = 0  no errors
    !           = 1  error in the specification of nlat
    !           = 2  error in the specification of nlon
    !           = 3  error in the specification of isym
    !           = 4  error in the specification of nt
    !           = 5  error in the specification of idv
    !           = 6  error in the specification of jdv
    !           = 7  error in the specification of mdb
    !           = 8  error in the specification of ndb
    !           = 9  error in the specification of lshsgc
    !
    module subroutine divgc(nlat, nlon, isym, nt, dv, idv, jdv, br, bi, mdb, ndb, &
        wshsgc, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: isym
        integer(ip), intent(in)  :: nt
        real(wp),    intent(out) :: dv(idv, jdv, nt)
        integer(ip), intent(in)  :: idv
        integer(ip), intent(in)  :: jdv
        real(wp),    intent(in)  :: br(mdb, ndb, nt)
        real(wp),    intent(in)  :: bi(mdb, ndb, nt)
        integer(ip), intent(in)  :: mdb
        integer(ip), intent(in)  :: ndb
        real(wp),    intent(in)  :: wshsgc(:)
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip) :: required_wavetable_size
        type(ScalarSynthesisUtility) :: util

        ! Check input arguments
        required_wavetable_size = util%get_lshsgc(nlat, nlon)

        call util%check_vector_transform_inputs(isym, idv, jdv, &
            mdb, ndb, nlat, nlon, nt, required_wavetable_size, &
            wshsgc, ierror)

        ! Check error flag
        if (ierror /= 0) return

        call divergence_lower_utility_routine(nlat, nlon, isym, nt, dv, &
            idv, jdv, br, bi, wshsgc, shsgc, ierror)

    end subroutine divgc

end submodule divergence_gaussian_grid
