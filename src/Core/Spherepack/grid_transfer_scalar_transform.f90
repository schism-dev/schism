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
submodule(grid_transfer_routines) grid_transfer_scalar_transform

contains
    !     subroutine trssph(intl, igrida, nlona, nlata, da, igridb, nlonb, nlatb, db, ierror)
    !
    ! Purpose
    !
    !     subroutine trssph transfers data given in array da on a grid on the
    !     full sphere to data in array db on a grid on the full sphere.  the
    !     grids on which da is given and db is generated can be specified
    !     independently of each other (see description below and the arguments
    !     igrida, igridb).  for transferring vector data on the sphere, use
    !     subroutine trvsph.

    !     notice that scalar and vector quantities are fundamentally different
    !     on the sphere.  for example, vectors are discontinuous and multiple
    !     valued at the poles.  scalars are continuous and single valued at the
    !     poles. erroneous results would be produced if one attempted to transfer
    !     vector fields between grids with subroutine trssph applied to each
    !     component of the vector.
    !
    !
    ! Underlying grid assumptions and a description
    !
    !     discussions with the ncar scd data support group and others indicate
    !     there is no standard grid for storing observational or model generated
    !     data on the sphere.  subroutine trssph was designed to handle most
    !     cases likely to be encountered when moving data from one grid format
    !     to another.
    !
    !     the grid on which da is given must be equally spaced in longitude
    !     and either equally spaced or gaussian in latitude (or colatitude).
    !     longitude, which can be either the first or second dimension of da,
    !     subdivides [0, 2pi) excluding the periodic point 2pi.  (co)latitude,
    !     which can be the second or first dimension of da, has south
    !     to north or north to south orientation with increasing subscript
    !     value in da (see the argument igrida).
    !
    !     the grid on which db is generated must be equally spaced in longitude
    !     and either equally spaced or gaussian in latitude (or colatitude).
    !     longitude, which can be either the first or second dimension of db,
    !     subdivides [0, 2pi) excluding the periodic point 2pi.  (co)latitude,
    !     which can be the second or first dimension of db, has south
    !     to north or north to south orientation with increasing subscript
    !     value in db (see the argument igridb).
    !
    !     let nlon be either nlona or nlonb (the number of grid points in
    !     longitude.  the longitude grid subdivides [0, 2pi) into nlon spaced
    !     points
    !
    !          (j-1)*2.*pi/nlon  (j=1, ..., nlon).
    !
    !     it is not necessary to communicate to subroutine trssph whether the
    !     underlying grids are in latitude or colatitude.  it is only necessary
    !     to communicate whether they run south to north or north to south with
    !     increasing subscripts.  a brief discussion of latitude and colatitude
    !     follows.  equally spaced latitude grids are assumed to subdivide
    !     [-pi/2, pi/2] with the south pole at -pi/2 and north pole at pi/2.
    !     equally spaced colatitude grids subdivide [0, pi] with the north pole
    !     at 0 and south pole at pi.  equally spaced partitions on the sphere
    !     include both poles.  gaussian latitude grids subdivide (-pi/2, pi/2)
    !     and gaussian colatitude grids subdivide (0, pi).  gaussian grids do not
    !     include the poles.  the gaussian grid points are uniquely determined by
    !     the size of the partition.  they can be computed in colatitude in
    !     (0, pi) (north to south) in real by the spherepack subroutine
    !     compute_gaussian_latitudes_and_weights.  let nlat be nlata or nlatb if either the da or db grid is
    !     gaussian.  let
    !
    !        north pole                             south pole
    !        ----------                             ----------
    !           0.0    <  cth(1) < ... < cth(nlat)  <   pi
    !
    !
    !     be nlat gaussian colatitude points in the interval (0, pi) and let
    !
    !        south pole                        north pole
    !        ----------                        ----------
    !           -pi/2  < th(1) < ... < th(nlat) < pi/2
    !
    !     be nlat gaussian latitude points in the open interval (-pi/2, pi/2).
    !     these are related by
    !
    !          th(i) = -pi/2 + cth(i)  (i=1, ..., nlat)
    !
    !     if the da or db grid is equally spaced in (co)latitude then
    !
    !          ctht(i) = (i-1)*pi/(nlat-1)
    !                                               (i=1, ..., nlat)
    !          tht(i) = -pi/2 + (i-1)*pi/(nlat-1)
    !
    !     define the equally spaced (north to south) colatitude and (south to
    !     north) latitude grids.
    !
    !
    ! Method (simplified description)
    !
    !     for simplicity, assume da is a nlat by nlon data tabulation and da(i, j)
    !     is the value at latitude theta(i) and longitude phi(j).  then
    !     coefficients a(m, n) and b(m, n) can be determined so that da(i, j) is
    !     approximated by the sum
    !
    !         l-1  n
    !     (a) sum sum pbar(m, n, theta(i))*(a(m, n)*cos(m*phi(j)+b(m, n)*sin(m*phi(j))
    !         n=0 m=0
    !
    !     here pbar(n, m, theta) are the normalized associated legendre functions
    !     and l = min(nlat, (nlon+2)/2).  the determination of a(m, n) and b(m, n)
    !     is called spherical harmonic analysis. a sum of this form can then be
    !     used to regenerate the data in db on the new grid with the known
    !     a(m, n) and b(m, n).  this is referred to spherical harmonic synthesis.
    !     analysis and synthesis subroutines from the software package spherepack,
    !     are used for these purposes.
    !
    !     if da or db is not in mathematical spherical coordinates then array
    !     transposition and/or subscript reordering is used prior to harmonic
    !     analysis and after harmonic synthesis.
    !
    ! Advantages
    !
    !     the use of surface spherical harmonics to transfer spherical grid data
    !     has advantages over pointwise grid interpolation schemes on the sphere.
    !     it is highly accurate.  if p(x, y, z) is any polynomial of degree n or
    !     less in x, y, z cartesian coordinates which is restricted to the surface
    !     of the sphere, then p is exactly represented by sums of the form (a)
    !     whenever n = mino(nlat, nlon/2) (i.e., transfers with spherical harmonics
    !     have n(th) order accuracy.  by way of contrast, bilinear interpolation
    !     schemes are exact for polynomials of degree one.  bicubic interpolation
    !     is exact only for polynomials of degree three or less.  the method
    !     also produces a weighted least squares fit to the data in which waves
    !     are resolved uniformly on the full sphere.  high frequencies, induced
    !     by closeness of grid points near the poles (due to computational
    !     or observational errors) are smoothed.  finally, the method is
    !     consistent with methods used to generate data in numerical spectral
    !     models based on spherical harmonics.  for more discussion of these and
    !     related issues,  see the article: "on the spectral approximation of
    !     discrete scalar and vector functions on the sphere, " siam j. numer.
    !     anal., vol 16. dec 1979, pp. 934-949, by paul swarztrauber.
    !
    !
    ! Comment
    !
    !     on a nlon by nlat or nlat by nlon grid (gaussian or equally spaced)
    !     spherical harmonic analysis generates and synthesis utilizes
    !     min(nlat, (nlon+2)/2)) by nlat coefficients.  consequently, for
    !     da and db,  if either
    !
    !             min(nlatb, (nlonb+2)/2) < min(nlata, (nlona+2)/2)
    !
    !     or if
    !
    !             nlatb < nlata
    !
    !     then all the coefficients generated by an analysis of da cannot be used
    !     in the synthesis which generates db.  in this case "information" can be
    !     lost in generating db.  more precisely, information will be lost if the
    !     analysis of da yields nonzero coefficients which are outside the bounds
    !     determined by the db grid.  nevertheless, transference of values with
    !     spherical harmonics will yield results consistent with grid resolution
    !     and is highly accurate.
    !
    !
    ! Input arguments
    !
    ! ... intl
    !
    !     an initialization argument which should be zero on an initial call to
    !     trssph.  intl should be one if trssph is being recalled and
    !
    !          igrida, nlona, nlata, igridb, nlonb, nlatb
    !
    !     have not changed from the previous call.  if any of these arguments
    !     have changed, intl=0 must be used to avoid undetectable errors.  calls
    !     with intl=1 bypass redundant computation and save time.  it can be used
    !     when transferring multiple data sets with the same underlying grids.
    !
    !
    ! ... igrida
    !
    !     an integer vector dimensioned two which identifies the underlying grid
    !     on the full sphere for the given data array da as follows:
    !
    !     igrida(1)
    !
    !     = -1
    !     if the latitude (or colatitude) grid for da is an equally spaced
    !     partition of [-pi/2, pi/2] ( or [0, pi]) including the poles which
    !     runs north to south
    !
    !     = +1
    !     if the latitude (or colatitude) grid for da is an equally spaced
    !     partition of [-pi/2, pi/2] ( or [0, pi]) including the poles which
    !     runs south to north
    !
    !     = -2
    !     if the latitude (or colatitude) grid for da is a gaussian partition
    !     of (-pi/2, pi/2) ( or (0, pi)) excluding the poles which runs north
    !     to south
    !
    !     = +2
    !     if the latitude (or colatitude) grid for da is a gaussian partition
    !     of (-pi/2, pi/2) ( or (0, pi)) excluding the poles which runs south
    !     north
    !
    !     igrida(2)
    !
    !     = 0 if the underlying grid for da is a nlona by nlata
    !
    !     = 1 if the underlying grid for da is a nlata by nlona
    !
    !
    ! ... nlona
    !
    !     the number of longitude points on the uniform grid which partitions
    !     [0, 2pi) for the given data array da.  nlona is also the first or second
    !     dimension of da (see igrida(2)) in the program which calls trssph.
    !     nlona determines the grid increment in longitude as 2*pi/nlona. for
    !     example nlona = 72 for a five degree grid.  nlona must be greater than
    !     or equal to 4.  the efficiency of the computation is improved when
    !     nlona is a product of small prime numbers
    !
    ! ... nlata
    !
    !     the number of points in the latitude (or colatitude) grid
    !     for the given data array da.  nlata is also the first or second
    !     dimension of da (see igrida(2)) in the program which calls trssph.
    !     if nlata is odd then the equator will be located at the (nlata+1)/2
    !     gaussian grid point.  if nlata is even then the equator will be
    !     located half way between the nlata/2 and nlata/2+1 grid points.
    !
    ! Note:
    !     igrida(1)=-1 or igrida(1)=-2 and igrida(2)=1 corresponds to
    !     the "usual" mathematical spherical coordinate system required
    !     by most of the drivers in spherepack2.  igrida(1)=1 or igrida(1)=2
    !     and igrida(2)=0 corresponds to the "usual" geophysical spherical
    !     coordinate system.
    !
    ! ... da
    !
    !     a two dimensional array that contains the data to be transferred.
    !     da must be dimensioned nlona by nlata in the program calling trssph if
    !     igrida(2) = 0.  da must be dimensioned nlata by nlona in the program
    !     calling trssph if igrida(2) = 1.  if da is not properly dimensioned
    !     and if the latitude (colatitude) values do not run south to north or
    !     north to south as flagged by igrida(1) (self cannot be checked!) then
    !     incorrect results will be produced.
    !
    ! ... igridb
    !
    !     an integer vector dimensioned two which identifies the underlying grid
    !     on the full sphere for the transformed data array db as follows:
    !
    !     igridb(1)
    !
    !     = -1
    !     if the latitude (or colatitude) grid for db is an equally spaced
    !     partition of [-pi/2, pi/2] ( or [0, pi]) including the poles which
    !     north to south
    !
    !     = +1
    !     if the latitude (or colatitude) grid for db is an equally spaced
    !     partition of [-pi/2, pi/2] ( or [0, pi]) including the poles which
    !     south to north
    !
    !     = -2
    !     if the latitude (or colatitude) grid for db is a gaussian partition
    !     of (-pi/2, pi/2) ( or (0, pi)) excluding the poles which runs north to
    !     south
    !
    !     = +2
    !     if the latitude (or colatitude) grid for db is a gaussian partition
    !     of (-pi/2, pi/2) ( or (0, pi)) excluding the poles which runs south to
    !     north
    !
    !
    !     igridb(2)
    !
    !     = 0 if the underlying grid for db is a nlonb by nlatb
    !
    !     = 1 if the underlying grid for db is a nlatb by nlonb
    !
    !
    ! ... nlonb
    !
    !     the number of longitude points on the uniform grid which partitions
    !     [0, 2pi) for the transformed data array db.  nlonb is also the first or
    !     second dimension of db (see igridb(2)) in the program which calls
    !     trssph.  nlonb determines the grid increment in longitude as 2*pi/nlonb.
    !     for example nlonb = 72 for a five degree grid.  nlonb must be greater
    !     than or equal to 4.  the efficiency of the computation is improved when
    !     nlonb is a product of small prime numbers
    !
    ! ... nlatb
    !
    !     the number of points in the latitude (or colatitude) grid
    !     for the transformed data array db.  nlatb is also the first or second
    !     dimension of db (see igridb(2)) in the program which calls trssph.
    !     if nlatb is odd then the equator will be located at the (nlatb+1)/2
    !     gaussian grid point.  if nlatb is even then the equator will be
    !     located half way between the nlatb/2 and nlatb/2+1 grid points.
    !
    ! Output arguments
    !
    !
    ! ... db
    !
    !     a two dimensional array that contains the transformed data.  db
    !     must be dimensioned nlonb by nlatb in the program calling trssph if
    !     igridb(2) = 0 or 1.  db must be dimensioned nlatb by nlonb in the
    !     program calling trssph if igridb(2) = 1.  if db is not properly
    !     dimensioned and if the latitude (colatitude) values do not run south
    !     north or north to south as flagged by igrdb(1) (self cannot be checked!)
    !     then incorrect results will be produced.
    !
    ! Error argument
    !
    ! ... ierror = 0  if no errors are detected
    !
    !         = 1  if intl is not 0 or 1
    !
    !         = 2  if igrida(1) is not -1 or +1 or -2 or +2
    !
    !         = 3  if igrida(2) is not 0 or 1
    !
    !         = 4  if nlona is less than 4
    !
    !         = 5  if nlata is less than 3
    !
    !         = 6  if igridb(1) is not -1 or +1 or -2 or +2
    !
    !         = 7  if igridb(2) is not 0 or 1
    !
    !         = 8  if nlonb is less than 4
    !
    !         = 9  if nlatb is less than 3
    !
    !         =12  indicates failure in an eigenvalue routine which computes
    !              gaussian weights and points
    !
    module subroutine trssph(intl, igrida, nlona, nlata, da, igridb, nlonb, nlatb, db, ierror)

        ! Dummy arguments
        integer(ip), intent(in)     :: intl
        integer(ip), intent(in)     :: igrida(2)
        integer(ip), intent(in)     :: nlona
        integer(ip), intent(in)     :: nlata
        real(wp),    intent(inout)  :: da(:,:)
        integer(ip), intent(in)     :: igridb(2)
        integer(ip), intent(in)     :: nlonb
        integer(ip), intent(in)     :: nlatb
        real(wp),    intent(out)    :: db(:,:)
        integer(ip), intent(out)    :: ierror

        ! Local variables
        integer(ip)                  :: igrda, igrdb, mdab_a, mdab_b
        real(wp), allocatable        :: wavetable_a(:), wavetable_b(:)
        integer(ip), parameter       :: NT = 1, ISYM = 0
        type(ScalarForwardTransform)  :: analysis_util
        type(ScalarBackwardTransform) :: synthesis_util

        ! Check calling arguments
        associate( &
            ig1 => igrida(1), &
            ig2 => igrida(2) &
            )
            if (intl*(intl-1)/=0) then
                ierror = 1
            else if ((ig1-1)*(ig1+1)*(ig1-2)*(ig1+2)/=0) then
                ierror = 2
            else if (ig2*(ig2-1)/=0) then
                ierror = 3
            else if (nlona < 4) then
                ierror = 4
            else if (nlata <3) then
                ierror = 5
            else if ((ig1-1)*(ig1+1)*(ig1-2)*(ig1+2) /= 0) then
                ierror = 6
            else if (ig2*(ig2-1) /= 0) then
                ierror = 7
            else if (nlonb < 4) then
                ierror = 8
            else if (nlatb < 3) then
                ierror = 9
            else
                ierror = 0
            end if
        end associate

        igrda = abs(igrida(1))
        igrdb = abs(igridb(1))

        if (igrda == 1) then
            ! Initialize wavetable for equally spaced analysis
            call analysis_util%initialize_shaec(nlata, nlona, wavetable_a, ierror)
        else
            ! Initialize wavetable for gaussian analysis
            call analysis_util%initialize_shagc(nlata, nlona, wavetable_a, ierror)
        end if

        if (igrdb == 2) then
            ! Initialize wavetable for gaussian synthesis
            call synthesis_util%initialize_shsgc(nlatb, nlonb, wavetable_b, ierror)
        else
            ! Initialize wsave for equally spaced synthesis
            call synthesis_util%initialize_shsec(nlatb, nlonb, wavetable_b, ierror)
        end if

        ! Transpose and/or reorder (co)latitude if necessary for da
        ! (arrays must have latitude (colatitude) as the first dimension
        ! and run north to south for spherepack software)
        if (igrida(2) == 0) call transpose_array(nlona, nlata, da)
        if (igrida(1) > 0) call reverse_colatitudes(nlata, nlona, da)

        mdab_a = min(nlata, (nlona + 2)/2)
        mdab_b = min(nlatb, (nlonb + 2)/2)

        block
            real(wp), dimension(mdab_a, nlata, NT) :: br_a, bi_a
            real(wp), dimension(mdab_b, nlatb, NT) :: br_b, bi_b

            if (igrda == 2) then
                ! Spherical harmonic analysis of "adjusted" da on gaussian grid
                call analysis_util%shagc(nlata, nlona, ISYM, NT, da, nlata, nlona, br_a, &
                    bi_a, mdab_a, nlata, wavetable_a, ierror)
            else
                ! Spherical harmonic analysis of "adjusted" da on equally spaced grid
                call analysis_util%shaec(nlata, nlona, ISYM, NT, da, nlata, nlona, br_a, &
                    bi_a, mdab_a, nlata, wavetable_a, ierror)
            end if

            ! Transfer da grid coefficients to db grid coefficients
            ! truncating to zero as necessary
            call transfer_scalar_coeff(mdab_a, nlata, br_a, bi_a, mdab_b, nlatb, br_b, &
                bi_b)

            if (igrdb == 1) then
                ! Spherical harmonic synthesis on nlatb by nlonb equally spaced grid
                call synthesis_util%shsec(nlatb, nlonb, ISYM, NT, db, nlatb, nlonb, br_b, &
                    bi_b, mdab_b, nlatb, wavetable_b, ierror)
            else
                ! Spherical harmonic synthesis on nlatb by nlonb gaussian grid
                call synthesis_util%shsgc(nlatb, nlonb, ISYM, NT, db, nlatb, nlonb, br_b, &
                    bi_b, mdab_b, nlatb, wavetable_b, ierror)
            end if
        end block

        ! Both da, db are currently latitude by longitude north to south arrays
        ! restore da and set db to agree with flags in igrida and igridb
        if (igrida(1) > 0) call reverse_colatitudes(nlata, nlona, da)

        if (igridb(1) > 0) call reverse_colatitudes(nlatb, nlonb, db)

        if (igrida(2) == 0) call transpose_array(nlata, nlona, da)

        if (igridb(2) == 0) call transpose_array(nlatb, nlonb, db)

        ! Release memory
        deallocate (wavetable_a)
        deallocate (wavetable_b)

    end subroutine trssph

end submodule grid_transfer_scalar_transform
