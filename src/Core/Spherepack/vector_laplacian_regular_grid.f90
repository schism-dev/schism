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
submodule(vector_laplacian_routines) vector_laplacian_regular_grid

contains
    !
    !
    ! ... file vlapec.f
    !
    !     this file includes documentation and code for
    !     subroutine vlapec          i
    !
    ! ... files which must be loaded with vlapec.f
    !
    !     type_SpherepackUtility.f, type_RealPeriodicFastFourierTransform.f, vhaec.f, vhsec.f
    !
    !
    !     subroutine vlapec(nlat, nlon, ityp, nt, vlap, wlap, idvw, jdvw, br, bi, cr, ci,
    !    +mdbc, ndbc, wvhsec, lvhsec, work, lwork, ierror)
    !
    !
    !     subroutine vlapec computes the vector laplacian of the vector field
    !     (v, w) in (vlap, wlap) (see the definition of the vector laplacian at
    !     the output parameter description of vlap, wlap below).  w and wlap
    !     are east longitudinal components of the vectors.  v and vlap are
    !     colatitudinal components of the vectors.  br, bi, cr, and ci are the
    !     vector harmonic coefficients of (v, w).  these must be precomputed by
    !     vhaec and are input parameters to vlapec.  the laplacian components
    !     in (vlap, wlap) have the same symmetry or lack of symmetry about the
    !     equator as (v, w).  the input parameters ityp, nt, mdbc, nbdc must have
    !     the same values used by vhaec to compute br, bi, cr, and ci for (v, w).
    !     vlap(i, j) and wlap(i, j) are given on the sphere at the colatitude
    !
    !            theta(i) = (i-1)*pi/(nlat-1)
    !
    !     for i=1, ..., nlat and east longitude
    !
    !            lambda(j) = (j-1)*2*pi/nlon
    !
    !     for j=1, ..., nlon.
    !
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
    !     ityp   this parameter should have the same value input to subroutine
    !            vhaec to compute the coefficients br, bi, cr, and ci for the
    !            vector field (v, w).  ityp is set as follows:
    !
    !            = 0  no symmetries exist in (v, w) about the equator. (vlap, wlap)
    !                 is computed and stored on the entire sphere in the arrays
    !                 vlap(i, j) and wlap(i, j) for i=1, ..., nlat and j=1, ..., nlon.
    !
    !
    !            = 1  no symmetries exist in (v, w) about the equator. (vlap, wlap)
    !                 is computed and stored on the entire sphere in the arrays
    !                 vlap(i, j) and wlap(i, j) for i=1, ..., nlat and j=1, ..., nlon.
    !                 the vorticity of (v, w) is zero so the coefficients cr and
    !                 ci are zero and are not used.  the vorticity of (vlap, wlap)
    !                 is also zero.
    !
    !
    !            = 2  no symmetries exist in (v, w) about the equator. (vlap, wlap)
    !                 is computed and stored on the entire sphere in the arrays
    !                 vlap(i, j) and wlap(i, j) for i=1, ..., nlat and j=1, ..., nlon.
    !                 the divergence of (v, w) is zero so the coefficients br and
    !                 bi are zero and are not used.  the divergence of (vlap, wlap)
    !                 is also zero.
    !
    !            = 3  w is antisymmetric and v is symmetric about the equator.
    !                 consequently wlap is antisymmetric and vlap is symmetric.
    !                 (vlap, wlap) is computed and stored on the northern
    !                 hemisphere only.  if nlat is odd, storage is in the arrays
    !                 vlap(i, j), wlap(i, j) for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon.
    !                 if nlat is even, storage is in the arrays vlap(i, j),
    !                 wlap(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !            = 4  w is antisymmetric and v is symmetric about the equator.
    !                 consequently wlap is antisymmetric and vlap is symmetric.
    !                 (vlap, wlap) is computed and stored on the northern
    !                 hemisphere only.  if nlat is odd, storage is in the arrays
    !                 vlap(i, j), wlap(i, j) for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon.
    !                 if nlat is even, storage is in the arrays vlap(i, j),
    !                 wlap(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.  the
    !                 vorticity of (v, w) is zero so the coefficients cr, ci are
    !                 zero and are not used. the vorticity of (vlap, wlap) is
    !                 also zero.
    !
    !            = 5  w is antisymmetric and v is symmetric about the equator.
    !                 consequently wlap is antisymmetric and vlap is symmetric.
    !                 (vlap, wlap) is computed and stored on the northern
    !                 hemisphere only.  if nlat is odd, storage is in the arrays
    !                 vlap(i, j), wlap(i, j) for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon.
    !                 if nlat is even, storage is in the arrays vlap(i, j),
    !                 wlap(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.  the
    !                 divergence of (v, w) is zero so the coefficients br, bi
    !                 are zero and are not used. the divergence of (vlap, wlap)
    !                 is also zero.
    !
    !
    !            = 6  w is symmetric and v is antisymmetric about the equator.
    !                 consequently wlap is symmetric and vlap is antisymmetric.
    !                 (vlap, wlap) is computed and stored on the northern
    !                 hemisphere only.  if nlat is odd, storage is in the arrays
    !                 vlap(i, j), wlap(i, j) for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon.
    !                 if nlat is even, storage is in the arrays vlap(i, j),
    !                 wlap(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !            = 7  w is symmetric and v is antisymmetric about the equator.
    !                 consequently wlap is symmetric and vlap is antisymmetric.
    !                 (vlap, wlap) is computed and stored on the northern
    !                 hemisphere only.  if nlat is odd, storage is in the arrays
    !                 vlap(i, j), wlap(i, j) for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon.
    !                 if nlat is even, storage is in the arrays vlap(i, j),
    !                 wlap(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.  the
    !                 vorticity of (v, w) is zero so the coefficients cr, ci are
    !                 zero and are not used. the vorticity of (vlap, wlap) is
    !                 also zero.
    !
    !            = 8  w is symmetric and v is antisymmetric about the equator.
    !                 consequently wlap is symmetric and vlap is antisymmetric.
    !                 (vlap, wlap) is computed and stored on the northern
    !                 hemisphere only.  if nlat is odd, storage is in the arrays
    !                 vlap(i, j), wlap(i, j) for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon.
    !                 if nlat is even, storage is in the arrays vlap(i, j),
    !                 wlap(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.  the
    !                 divergence of (v, w) is zero so the coefficients br, bi
    !                 are zero and are not used. the divergence of (vlap, wlap)
    !                 is also zero.
    !
    !
    !     nt     nt is the number of vector fields (v, w).  some computational
    !            efficiency is obtained for multiple fields.  in the program
    !            that calls vlapec, the arrays vlap, wlap, br, bi, cr and ci
    !            can be three dimensional corresponding to an indexed multiple
    !            vector field.  in this case multiple vector synthesis will
    !            be performed to compute the vector laplacian for each field.
    !            the third index is the synthesis index which assumes the values
    !            k=1, ..., nt.  for a single synthesis set nt=1.  the description
    !            of the remaining parameters is simplified by assuming that nt=1
    !            or that all arrays are two dimensional.
    !
    !   idvw     the first dimension of the arrays vlap and wlap as it appears
    !            in the program that calls vlapec.  if ityp=0, 1, or 2  then idvw
    !            must be at least nlat.  if ityp > 2 and nlat is even then idvw
    !            must be at least nlat/2. if ityp > 2 and nlat is odd then idvw
    !            must be at least (nlat + 1)/2.
    !
    !   jdvw     the second dimension of the arrays vlap and wlap as it appears
    !            in the program that calls vlapec. jdvw must be at least nlon.
    !
    !
    !   br, bi    two or three dimensional arrays (see input parameter nt)
    !   cr, ci    that contain vector spherical harmonic coefficients
    !            of the vector field (v, w) as computed by subroutine vhaec.
    !            br, bi, cr and ci must be computed by vhaec prior to calling
    !            vlapec.  if ityp=1, 4, or 7 then cr, ci are not used and can
    !            be dummy arguments.  if ityp=2, 5, or 8 then br, bi are not
    !            used and can be dummy arguments.
    !
    !    mdbc    the first dimension of the arrays br, bi, cr and ci as it
    !            appears in the program that calls vlapec.  mdbc must be
    !            at least min(nlat, nlon/2) if nlon is even or at least
    !            min(nlat, (nlon + 1)/2) if nlon is odd.
    !
    !    ndbc    the second dimension of the arrays br, bi, cr and ci as it
    !            appears in the program that calls vlapec. ndbc must be at
    !            least nlat.
    !
    !    wvhsec  an array which must be initialized by subroutine vhseci.
    !            once initialized, wvhsec
    !            can be used repeatedly by vlapec as long as nlat and nlon
    !            remain unchanged.  wvhsec must not be altered between calls
    !            of vlapec.
    !
    !    lvhsec  the dimension of the array wvhsec as it appears in the
    !            program that calls vlapec.  let
    !
    !               l1 = min(nlat, (nlon+2)/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd.
    !
    !            then lvhsec must be at least
    !
    !            4*nlat*l2+3*max(l1-2, 0)*(2*nlat-l1-1)+nlon+15
    !
    !     output parameters
    !
    !
    !    vlap,   two or three dimensional arrays (see input parameter nt) that
    !    wlap    contain the vector laplacian of the field (v, w).  wlap(i, j) is
    !            the east longitude component and vlap(i, j) is the colatitudinal
    !            component of the vector laplacian.  the definition of the
    !            vector laplacian follows:
    !
    !            let cost and sint be the cosine and sine at colatitude theta.
    !            let d()/dlambda  and d()/dtheta be the first order partial
    !            derivatives in longitude and colatitude.  let del2 be the scalar
    !            laplacian operator
    !
    !                 del2(s) = [d(sint*d(s)/dtheta)/dtheta +
    !                             2            2
    !                            d (s)/dlambda /sint]/sint
    !
    !            then the vector laplacian opeator
    !
    !                 dvel2(v, w) = (vlap, wlap)
    !
    !            is defined by
    !
    !                 vlap = del2(v) - (2*cost*dw/dlambda + v)/sint**2
    !
    !                 wlap = del2(w) + (2*cost*dv/dlambda - w)/sint**2
    !
    !  ierror    a parameter which flags errors in input parameters as follows:
    !
    !            = 0  no errors detected
    !            = 1  error in the specification of nlat
    !            = 2  error in the specification of nlon
    !            = 3  error in the specification of ityp
    !            = 4  error in the specification of nt
    !            = 5  error in the specification of idvw
    !            = 6  error in the specification of jdvw
    !            = 7  error in the specification of mdbc
    !            = 8  error in the specification of ndbc
    !            = 9  error in the specification of lvhsec
    !
    module subroutine vlapec(nlat, nlon, ityp, nt, vlap, wlap, idvw, jdvw, br, bi, &
        cr, ci, mdbc, ndbc, wvhsec, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: ityp
        integer(ip), intent(in)  :: nt
        real(wp),    intent(out) :: vlap(idvw, jdvw, nt)
        real(wp),    intent(out) :: wlap(idvw, jdvw, nt)
        integer(ip), intent(in)  :: idvw
        integer(ip), intent(in)  :: jdvw
        real(wp),    intent(in)  :: br(mdbc, ndbc, nt)
        real(wp),    intent(in)  :: bi(mdbc, ndbc, nt)
        real(wp),    intent(in)  :: cr(mdbc, ndbc, nt)
        real(wp),    intent(in)  :: ci(mdbc, ndbc, nt)
        integer(ip), intent(in)  :: mdbc
        integer(ip), intent(in)  :: ndbc
        real(wp),    intent(in)  :: wvhsec(:)
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip) :: required_wavetable_size
        type(VectorSynthesisUtility) :: util
        
        ! Check input arguments
        required_wavetable_size = util%get_lvhsec(nlat, nlon)

        call util%check_vector_transform_inputs(ityp, idvw, jdvw, &
            mdbc, ndbc, nlat, nlon, nt, required_wavetable_size, &
            wvhsec, ierror)

        ! Check error flag
        if (ierror /= 0) return

        call vector_laplacian_lower_utility_routine(nlat, nlon, ityp, nt, vlap, wlap, &
            br, bi, cr, ci, wvhsec, vhsec, ierror)

    end subroutine vlapec

end submodule vector_laplacian_regular_grid
