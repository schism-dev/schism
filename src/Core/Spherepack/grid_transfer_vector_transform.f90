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
submodule(grid_transfer_routines) grid_transfer_vector_transform

contains

    !     subroutine trvsph(intl, igrida, nlona, nlata, iveca, ua, va, &
    !                       igridb, nlonb, nlatb, ivecb, ub, vb, ierror)
    !
    ! *** author
    !
    !     John C. Adams (NCAR 1997), email: johnad@ncar.ucar.edu
    !
    ! *** purpose
    !
    !     subroutine trvsph transfers vector data given in (ua, va) on a grid on
    !     the full sphere to vector data in (ub, vb) on a grid on the full sphere.
    !     the grids on which (ua, va) is given and (ub, vb) is generated can be
    !     specified independently of each other (see the input arguments igrida,
    !     igridb, iveca, ivecb).  ua and ub are the east longitudinal components of
    !     the given and transformed vector fields.  va is either the latitudinal
    !     or colatitudinal component of the given vector field (see iveca).
    !     vb is either the latitudinal or colatitudinal component of the
    !     transformed vector field (see ivecb).  for transferring scalar data
    !     on the sphere, use subroutine trssph.
    !
    ! *   notice that scalar and vector quantities are fundamentally different
    !     on the sphere.  for example, vectors are discontinuous and multiple
    !     valued at the poles.  scalars are continuous and single valued at the
    !     poles. erroneous results would be produced if one attempted to transfer
    !     vector fields between grids with subroutine trssph applied to each
    !     component of the vector.
    !
    ! *** underlying grid assumptions and a description
    !
    !     discussions with the ncar scd data support group and others indicate
    !     there is no standard grid for storing observational or model generated
    !     data on the sphere.  subroutine trvsph was designed to handle most
    ! cases likely to be encountered when moving data from one grid format
    !     to another.
    !
    !     the grid on which (ua, va) is given must be equally spaced in longitude
    !     and either equally spaced or gaussian in latitude (or colatitude).
    !     longitude, which can be either the first or second dimension of ua, va
    !     subdivides [0, 2pi) excluding the periodic point 2pi.  (co)latitude,
    !     which can be the second or first dimension of ua, va, has south
    !     to north or north to south orientation with increasing subscript
    !     value in ua, va (see the argument igrida).
    !
    !     the grid on which ub, vb is generated must be equally spaced in longitude
    !     and either equally spaced or gaussian in latitude (or colatitude).
    !     longitude, which can be either the first or second dimension of ub, vb
    !     subdivides [0, 2pi) excluding the periodic point 2pi.  (co)latitude,
    !     which can be the second or first dimension of ub, vb, has south
    !     to north or north to south orientation with increasing subscript
    !     value in db (see the argument igridb).
    !
    !     let nlon be either nlona or nlonb (the number of grid points in
    !     longitude.  the longitude grid subdivides [0, 2pi) into nlon spaced
    !     points
    !
    !          (j-1)*2.*pi/nlon  (j=1, ..., nlon).
    !
    !     it is not necessary to communicate to subroutine trvsph whether the
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
    !     compute_gaussian_latitudes_and_weights.  let nlat be nlata or nlatb if either the ua, va or ub, vb grid is
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
    !     if the (ua, va) or (ub, vb) grid is equally spaced in (co)latitude then
    !
    !          ctht(i) = (i-1)*pi/(nlat-1)
    !                                               (i=1, ..., nlat)
    !          tht(i) = -pi/2 + (i-1)*pi/(nlat-1)
    !
    !     define the equally spaced (north to south) colatitude and (south to
    !     north) latitude grids.
    !
    ! *** method (simplified description)
    !
    !    (1)
    !
    !     the vector field (ua, va) is reformated to a vector field in mathematical
    !     spherical coordinates using array transpositions, subscript reordering
    !     and negation of va as necessary (see arguments igrida, iveca).
    !
    !     (2)
    !
    !     a vector harmonic analysis is performed on the result from (1)
    !
    !     (3)
    !
    !     a vector harmonic synthesis is performed on the (ub, vb) grid
    !     using as many coefficients from (2) as possible (i.e., as
    !     as is consistent with the size of the ub, vb grid).
    !
    !     (4)
    !
    !     the vector field generated in (3) is transformed from mathematical
    !     spherical coordinates to the form flagged by ivecb and igridb in
    !     (ub, vb) using array transpositions, subscript reordering and negation
    !     as necessary
    !
    !
    ! *** advantages
    !
    !     the use of vector spherical harmonics to transfer vector data is
    !     highly accurate and preserves properties of vectors on the sphere.
    !     the method produces a weighted least squares fit to vector data in
    !     which waves are resolved uniformly on the full sphere.  high frequencies
    !     induced by closeness of grid points near the poles (due to computational
    !     or observational errors) are smoothed.  the method is consistent with
    !     methods used to generate vector data in numerical spectral models based
    !     on spherical harmonics.  for more discussion of these and related issues,
    !     see "on the spectral approximation of discrete scalar and vector
    !     functions on the sphere, " siam j. numer. anal., vol. 16, december 1979, 
    !     pp. 934-949, by paul swarztrauber.
    !
    !
    ! *** comment
    !
    !     on a nlon by nlat or nlat by nlon grid (gaussian or equally spaced)
    !     spherical harmonic analysis generates and synthesis utilizes
    !     min(nlat, (nlon+2)/2)) by nlat coefficients.  consequently, for
    !     ua, va and ub, vb,  if either
    !
    !             min(nlatb, (nlonb+2)/2) < min(nlata, (nlona+2)/2)
    !
    !     or if
    !
    !             nlatb < nlata
    !
    !     then all the coefficients generated by an analysis of ua, va cannot be
    !     used in the synthesis which generates ub, vb.  in this case "information"
    !     can be lost in generating ub, vb.  more precisely, information will be
    !     lost if the analysis of ua, va yields nonzero coefficients which are
    !     outside the coefficient bounds determined by the ub, vb grid. still
    !     transference with vector spherical harmonics will yield results
    !     consistent with grid resolution and is highly accurate.
    !
    ! *** input arguments
    !
    ! ... intl
    !
    !     an initialization argument which should be zero on an initial call to
    !     trvsph.  intl should be one if trvsph is being recalled and
    !
    !          igrida, nlona, nlata, iveca, igridb, nlonb, nlatb, ivecb
    !
    !     have not changed from the previous call.  if any of these arguments have
    !     changed intl=0 must be used to avoid undetectable errors.  when allowed,
    !     calls with intl=1 bypass redundant computation and save time.  it can
    !     be used when transferring multiple vector data sets with the same
    !     underlying grids.
    !
    ! ... igrida
    !
    !     an integer vector dimensioned two which identifies the underlying grid
    !     on the full sphere for the given vector data (ua, va) as follows:
    !
    !     igrida(1)
    !
    !     = -1
    !     if the latitude (or colatitude) grid for ua, va is an equally spaced
    !     partition of [-pi/2, pi/2] ( or [0, pi]) including the poles which
    !     runs north to south with increasing subscript value
    !
    !     = +1
    !     if the latitude (or colatitude) grid for ua, va is an equally spaced
    !     partition of [-pi/2, pi/2] ( or [0, pi]) including the poles which
    !     runs south to north with increasing subscript value
    !
    !     = -2
    !     if the latitude (or colatitude) grid for ua, va is a gaussian partition
    !     of (-pi/2, pi/2) ( or (0, pi)) excluding the poles which runs north
    !     to south with increasing subscript value
    !
    !     = +2
    !     if the latitude (or colatitude) grid for ua, va is a gaussian partition
    !     of (-pi/2, pi/2) ( or (0, pi)) excluding the poles which runs south
    !     north with increasing subscript value
    !
    !     igrida(2)
    !
    !     = 0 if the underlying grid for ua, va is a nlona by nlata
    !
    !     = 1 if the underlying grid for ua, va is a nlata by nlona
    !
    !
    ! ... nlona
    !
    !     the number of longitude points on the uniform grid which partitions
    !     [0, 2pi) for the given vector (ua, va).  nlona is also the first or second
    !     dimension of ua, va (see igrida(2)) in the program which calls trvsph.
    !     nlona determines the grid increment in longitude as 2*pi/nlona. for
    !     example nlona = 72 for a five degree grid.  nlona must be greater than
    !     or equal to 4.  the efficiency of the computation is improved when
    !     nlona is a product of small prime numbers
    !
    ! ... nlata
    !
    !     the number of points in the latitude (or colatitude) grid for the
    !     given vector (ua, va).  nlata is also the first or second dimension
    !     of ua and va (see igrida(2)) in the program which calls trvsph.
    !     if nlata is odd then the equator will be located at the (nlata+1)/2
    !     gaussian grid point.  if nlata is even then the equator will be
    !     located half way between the nlata/2 and nlata/2+1 grid points.
    !
    ! ... iveca
    !
    !     if iveca=0 is input then va is the latitudinal component of the
    !     given vector field. if iveca=1 then va is the colatitudinal
    !     compoenent of the given vector field.  in either case, ua must
    !     be the east longitudinal component of the given vector field.
    !
    ! *** note:
    !     igrida(1)=-1 or igrida(1)=-2, igrida(2)=1, and iveca=1 corresponds
    !     to the "usual" mathematical spherical coordinate system required
    !     by most of the drivers in spherepack2.  igrida(1)=1 or igrida(1)=2,
    !     igrida(2)=0, and iveca=0 corresponds to the "usual" geophysical
    !     spherical coordinate system.
    !
    !
    ! ... ua
    !
    !     ua is the east longitudinal component of the given vector field.
    !     ua must be dimensioned nlona by nlata in the program calling trvsph if
    !     igrida(2) = 0.  ua must be dimensioned nlata by nlona in the program
    !     calling trvsph if igrida(2) = 1.  if ua is not properly dimensioned
    !     and if the latitude (colatitude) values do not run south to north or
    !     north to south as flagged by igrida(1) (self cannot be checked!) then
    !     incorrect results will be produced.
    !
    !
    ! ... va
    !
    !     va is either the latitudinal or colatitudinal componenet of the
    !     given vector field (see iveca).  va must be dimensioned nlona by
    !     nlata in the program calling trvsph if igrida(2)=0.  va must be
    !     dimensioned nlata by nlona in the program calling trvsph if
    !     igrida(2)=1.  if va is not properly dimensioned or if the latitude
    !     (colatitude) values do not run south to north or north to south
    !     as flagged by igrida(1) (self cannot be checked!) then incorrect
    !     results will be produced.
    !
    ! ... igridb
    !
    !     an integer vector dimensioned two which identifies the underlying grid
    !     on the full sphere for the transformed vector (ub, vb) as follows:
    !
    !     igridb(1)
    !
    !     = -1
    !     if the latitude (or colatitude) grid for ub, vb is an equally spaced
    !     partition of [-pi/2, pi/2] ( or [0, pi]) including the poles which
    !     north to south
    !
    !     = +1
    !     if the latitude (or colatitude) grid for ub, vb is an equally spaced
    !     partition of [-pi/2, pi/2] ( or [0, pi]) including the poles which
    !     south to north
    !
    !     = -2
    !     if the latitude (or colatitude) grid for ub, vb is a gaussian partition
    !     of (-pi/2, pi/2) ( or (0, pi)) excluding the poles which runs north to
    !     south
    !
    !     = +2
    !     if the latitude (or colatitude) grid for ub, vb is a gaussian partition
    !     of (-pi/2, pi/2) ( or (0, pi)) excluding the poles which runs south to
    !     north
    !
    !     igridb(2)
    !
    !     = 0 if the underlying grid for ub, vb is a nlonb by nlatb
    !
    !     = 1 if the underlying grid for ub, vb is a nlatb by nlonb
    !
    !
    ! ... nlonb
    !
    !     the number of longitude points on the uniform grid which partitions
    !     [0, 2pi) for the transformed vector (ub, vb).  nlonb is also the first or
    !     second dimension of ub and vb (see igridb(2)) in the program which calls
    !     trvsph.  nlonb determines the grid increment in longitude as 2*pi/nlonb.
    !     for example nlonb = 72 for a five degree grid.  nlonb must be greater
    !     than or equal to 4.  the efficiency of the computation is improved when
    !     nlonb is a product of small prime numbers
    !
    ! ... nlatb
    !
    !     the number of points in the latitude (or colatitude) grid for the
    !     transformed vector (ub, vb).  nlatb is also the first or second dimension
    !     of ub and vb (see igridb(2)) in the program which calls trvsph.
    !     if nlatb is odd then the equator will be located at the (nlatb+1)/2
    !     gaussian grid point.  if nlatb is even then the equator will be
    !     located half way between the nlatb/2 and nlatb/2+1 grid points.
    !
    ! ... ivecb
    !
    !     if ivecb=0 is input then vb is the latitudinal component of the
    !     given vector field. if ivecb=1 then vb is the colatitudinal
    !     compoenent of the given vector field.  in either case, ub must
    !     be the east longitudinal component of the given vector field.
    !
    ! *** note:
    !     igridb(1)=-1 or igridb(1)=-2, igridb(2)=1, and ivecb=1 corresponds
    !     to the "usual" mathematical spherical coordinate system required
    !     by most of the drivers in spherepack2.  igridb(1)=1 or igridb(1)=2,
    !     igridb(2)=0, and ivecb=0 corresponds to the "usual" geophysical
    !     spherical coordinate system.
    !
    ! ... wsave
    !
    !     a saved workspace array that can be utilized repeatedly by trvsph
    !     as long as the arguments nlata, nlona, nlatb, nlonb remain unchanged.
    !     wsave is set by a intl=0 call to trvsph.  wsave must not be altered
    !     when trvsph is being recalled with intl=1.
    !
    ! ... lsave
    !
    !     the dimension of the workspace wsave as it appears in the program
    !     that calls trvsph.  the minimum required value of lsave for the
    !     current set of input arguments is set in the output argument lsvmin.
    !     it can be determined by calling trvsph with lsave=0 and printing lsvmin.
    !
    !          la1 = min(nlata, (nlona+1)/2), la2 = (nlata+1)/2
    !
    !          lb1 = min(nlatb, (nlonb+1)/2), lb2 = (nlatb+1)/2
    !
    !          lwa = 4*nlata*la2+3*max(la1-2, 0)*(2*nlata-la1-1)+la2+nlona+15
    !
    !          lwb = 4*nlatb*lb2+3*max(lb1-2, 0)*(2*nlatb-lb1-1)+nlonb+15
    !
    !      then
    !
    !          lsvmin = lwa + lwb
    !
    !      is the minimal required workspace length of wsave
    !
    !
    ! ... work
    !
    !     a work array that does not have to be preserved
    !
    ! ... lwork
    !
    !     the dimension of the array work as it appears in the program that
    !     calls trvsph. the minimum required value of lwork for the current
    !     set of input arguments is set in the output argument lwkmin.
    !     it can be determined by calling trvsph with lwork=0 and printing
    !     lwkmin.  an estimate for lwork follows.  let nlat = max(nlata, nlatb),
    !     nlon = max(nlona, nlonb) and l1 = min(nlat, (nlon+2)/2).  with these
    !     these definitions, the quantity
    !
    !            2*nlat*(8*l1 + 4*nlon + 3)
    !
    !     will suffice as a length for the unsaved workspace.  this formula
    !     may overestimate the required minimum value for lwork.  the exact
    !     minimum value can be predetermined by calling trvsph wtih lwork=0
    !     and printout of lwkmin.
    !
    ! ... dwork
    !
    !     a real work array that does not have to be preserved.
    !
    ! ... ldwork
    !
    !     the length of dwork in the routine calling trvsph
    !     Let
    !
    !       nlat = max(nlata, nlatb)
    !
    !     ldwork must be at least 2*nlat*(nlat + 1)+1
    !
    !
    ! *** output arguments
    !
    !
    ! ... ub
    !
    !     a two dimensional array that contains the east longitudinal component
    !     of the transformed vector data.  ub
    !     must be dimensioned nlonb by nlatb in the program calling trvsph if
    !     igridb(2)=0.  ub must be dimensioned nlatb by nlonb in the program
    !     calling trvsph if igridb(2)=1.  if ub is not properly dimensioned
    !     and if the latitude (colatitude) values do not run south to north or
    !     north to south as flagged by igrdb(1) (self cannot be checked!) then
    !     incorrect results will be produced.
    !
    !
    ! ... vb
    !
    !     a two dimensional array that contains the latitudinal or colatitudinal
    !     component of the transformed vector data (see ivecb).
    !     vb must be dimensioned nlonb by nlatb in the program calling trvsph if
    !     igridb(2)=0.  vb must be dimensioned nlatb by nlonb in the program
    !     calling trvsph if igridb(2)=1.  if vb is not properly dimensioned
    !     and if the latitude (colatitude) values do not run south to north or
    !     north to south as flagged by igrdb(1) (self cannot be checked!) then
    !     incorrect results will be produced.
    !
    ! ... lsvmin
    !
    !     the minimum length of the saved workspace in wsave.
    !     lsvmin is computed even if lsave < lsvmin (ierror = 10).
    !
    ! ... lwkmin
    !
    !     the minimum length of the unsaved workspace in work.
    !     lwkmin is computed even if lwork < lwkmin (ierror = 11).
    !
    !
    ! *** error argument
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
    !         = 6  if iveca is not 0 or 1
    !
    !         = 7  if igridb(1) is not -1 or +1 or -2 or +2
    !
    !         = 8  if igridb(2) is not 0 or 1
    !
    !         = 9  if nlonb is less than 4
    !
    !         = 10  if nlatb is less than 3
    !
    !         = 11  if ivecb is not 0 or 1
    !
    !         = 14  indicates failure in an eigenvalue routine which computes
    !              gaussian weights and points
    !
    module subroutine trvsph(intl, igrida, nlona, nlata, iveca, ua, va, &
        igridb, nlonb, nlatb, ivecb, ub, vb, ierror)

        integer(ip), intent(in)    :: intl
        integer(ip), intent(in)    :: igrida(2)
        integer(ip), intent(in)    :: nlona, nlata, iveca
        real(wp),    intent(inout) :: ua(:,:), va(:,:)
        integer(ip), intent(in)    :: igridb(2), nlonb, nlatb
        integer(ip), intent(in)    :: ivecb
        real(wp),    intent(out)   :: ub(:,:), vb(:,:)
        integer(ip), intent(out)   :: ierror

        ! Local variables
        integer(ip)                  :: igrda, igrdb, mdab_a, mdab_b
        real(wp), allocatable        :: wavetable_a(:), wavetable_b(:)
        integer(ip), parameter       :: NT = 1, ITYP = 0
        type(VectorForwardTransform)  :: analysis_util
        type(VectorBackwardTransform) :: synthesis_util

        ! Check calling arguments
        associate( &
            ig1 => igrida(1), &
            ig2 => igrida(2) &
            )

            if (intl*(intl-1)/=0) then
                ierror = 1
            else if ((ig1-1)*(ig1+1)*(ig1-2)*(ig1+2) /= 0) then
                ierror = 2
            else if (ig2*(ig2-1) /= 0) then
                ierror = 3
            else if (nlona < 4) then
                ierror = 4
            else if (nlata < 3) then
                ierror = 5
            else if (iveca*(iveca-1)/=0) then
                ierror = 6
            else if ((ig1-1)*(ig1+1)*(ig1-2)*(ig1+2)/=0) then
                ierror = 7
            else if (ig2*(ig2-1) /= 0) then
                ierror = 8
            else if (nlonb < 4) then
                ierror = 9
            else if (nlatb < 3) then
                ierror = 10
            else if (ivecb*(ivecb-1) /= 0) then
                ierror = 11
            else
                ierror = 0
            end if
        end associate

        ! Check error flag
        if (ierror /= 0) return

        igrda = abs(igrida(1))
        igrdb = abs(igridb(1))

        if (igrda == 1) then
            ! Initialize wavetable for equally spaced analysis
            call analysis_util%initialize_vhaec(nlata, nlona, wavetable_a, ierror)
        else
            ! Initialize wavetable for gaussian analysis
            call analysis_util%initialize_vhagc(nlata, nlona, wavetable_a, ierror)
        end if

        if (igrdb == 2) then
            ! Initialize wavetable for gaussian synthesis
            call synthesis_util%initialize_vhsgc(nlatb, nlonb, wavetable_b, ierror)
        else
            ! Initialize wavetable for equally spaced synthesis
            call synthesis_util%initialize_vhsec(nlatb, nlonb, wavetable_b, ierror)
        end if

        ! Address error flag
        if (ierror /= 0) then
            ierror = 14
            return
        end if

        ! Convert the vector field (ua, va) to mathematical spherical coordinates
        if (igrida(2) == 0) then
            call transpose_array(nlona, nlata, ua)
            call transpose_array(nlona, nlata, va)
        end if

        if (igrida(1) > 0) then
            call reverse_colatitudes(nlata, nlona, ua)
            call reverse_colatitudes(nlata, nlona, va)
        end if

        if (iveca == 0) va = -va !call negate_colatitudinal_vector_component(nlata, nlona, va)

        ! Analyze vector field
        mdab_a = min(nlata, (nlona + 1)/2)
        mdab_b = min(nlatb, (nlonb + 1)/2)

        block
            real(wp), dimension(mdab_a, nlata, NT) :: br_a, bi_a, cr_a, ci_a
            real(wp), dimension(mdab_b, nlatb, NT) :: br_b, bi_b, cr_b, ci_b

            if (igrda == 2) then
                call analysis_util%vhagc(nlata, nlona, ITYP, NT, va, ua, nlata, nlona, &
                    br_a, bi_a, cr_a, ci_a, mdab_a, nlata, wavetable_a, ierror)
            else
                call analysis_util%vhaec(nlata, nlona, ITYP, NT, va, ua, nlata, nlona, &
                    br_a, bi_a, cr_a, ci_a, mdab_a, nlata, wavetable_a, ierror)
            end if

            ! Transfer a grid coefficients to b grid coefficients
            call transfer_vector_coeff(mdab_a, nlata, br_a, bi_a, cr_a, ci_a, &
                mdab_b, nlatb, br_b, bi_b, cr_b, ci_b)

            ! Synthesize on b grid
            if (igrdb == 1) then
                call synthesis_util%vhsec(nlatb, nlonb, ITYP, NT, vb, ub, nlatb, nlonb, br_b, &
                    bi_b, cr_b, ci_b, mdab_b, nlatb, wavetable_b, ierror)
            else
                call synthesis_util%vhsgc(nlatb, nlonb, ITYP, NT, vb, ub, nlatb, nlonb, br_b, &
                    bi_b, cr_b, ci_b, mdab_b, nlatb, wavetable_b, ierror)
            end if
        end block

        ! Restore a grid and b grid vector fields (now in math coordinates) to
        !  agree with grid flags in igrida, iveca, igridb, ivecb
        if (iveca == 0) va = -va
        if (ivecb == 0) vb = -vb

        if (igrida(1) > 0) then
            call reverse_colatitudes(nlata, nlona, ua)
            call reverse_colatitudes(nlata, nlona, va)
        end if

        if (igridb(1) > 0) then
            call reverse_colatitudes(nlatb, nlonb, ub)
            call reverse_colatitudes(nlatb, nlonb, vb)
        end if

        if (igrida(2) == 0) then
            call transpose_array(nlata, nlona, ua)
            call transpose_array(nlata, nlona, va)
        end if

        if (igridb(2) == 0) then
            call transpose_array(nlatb, nlonb, ub)
            call transpose_array(nlatb, nlonb, vb)
        end if

        ! Release memory
        deallocate (wavetable_a)
        deallocate (wavetable_b)

    end subroutine trvsph

    subroutine transfer_vector_coeff(ma, na, abr, abi, acr, aci, mb, nb, bbr, bbi, bcr, bci)

        ! Dummy arguments
        integer(ip), intent(in)  :: ma, na, mb, nb
        real(wp),    intent(in)  :: abr(ma, na), abi(ma, na), acr(ma, na), aci(ma, na)
        real(wp),    intent(out) :: bbr(mb, nb), bbi(mb, nb), bcr(mb, nb), bci(mb, nb)

        call transfer_scalar_coeff(ma, na, abr, abi, mb, nb, bbr, bbi)
        call transfer_scalar_coeff(ma, na, acr, aci, mb, nb, bcr, bci)

    end subroutine transfer_vector_coeff

end submodule grid_transfer_vector_transform
