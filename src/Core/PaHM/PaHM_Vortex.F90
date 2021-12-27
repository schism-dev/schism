!----------------------------------------------------------------
!               M O D U L E   V O R T E X
!----------------------------------------------------------------
!> @file vortex.F90
!>
!> @brief
!>   
!>
!> @details
!>   
!>
!> @author Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>
!> @note Adopted from the ADCIRC source code.
!----------------------------------------------------------------

!Routines & functions
! CalcIntensityChange
! UVTrans
! UVTransPoint
! NewVortex
! NewVortexFull
! SetVortex
! SetRMaxes
! GetRMaxes
! CalcRMaxes
! CalcRMaxesFull
! FitRMaxes
! FitRMaxes4
! SetVMaxesBL
! GetVMaxesBL
! SetUseVMaxesBL
! SetShapeParameter
! GetShapeParameter (fun)
! GetShapeParameters (fun)
! GetPhiFactors (fun)
! SetIsotachRadii
! SetIsotachWindSpeeds
! SetUsequadrantVR
! GetUsequadrantVR (fun)
! SpInterp (fun)
! InterpR (fun)
! Rmw (fun)
! UVP
! UVPR
! FAng (fun)
! Rotate
! GetLatestRMax (fun)
! GetLatestAngle (fun)
! VhWithCoriFull (fun)
! VhWithCori (fun)
! VhNoCori (fun)
! FindRoot (fun)


MODULE PaHM_Vortex

  USE PaHM_Sizes
  USE PaHM_Messages

  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER              :: NQUADS  = 4            ! Number of quadrants for which wind radii are provided
  INTEGER, PARAMETER              :: NPOINTS = NQUADS + 2   ! Number of (theta, rMax) points for curve fit
  REAL(SZ), DIMENSION(NPOINTS)    :: rMaxes                 ! Radius of maximum winds
  REAL(SZ), DIMENSION(NPOINTS, 4) :: rMaxes4                ! (nautical miles)

  REAL(SZ)                        :: pn                     ! Ambient surface pressure (mb) !PV global var?
  REAL(SZ)                        :: pc                     ! Surface pressure at center of storm (mb) !PV global var?
  REAL(SZ)                        :: cLat                   ! Latitude  of storm center (degrees north) !PV global var?
  REAL(SZ)                        :: cLon                   ! Longitude of storm center (degrees east ) !PV global var?
  REAL(SZ)                        :: vMax                   ! Max sustained wind velocity in storm (knots) !PV global var?

  REAL(SZ)                        :: B                      ! Exponential shape parameter
  REAL(SZ)                        :: corio                  ! Coriolis force (1/s)
  REAL(SZ)                        :: vr                     ! Velocity @ wind radii (knots)
  REAL(SZ)                        :: phi
  REAL(SZ), DIMENSION(NPOINTS)    :: phis                   ! Correction factor to B and vh
  REAL(SZ), DIMENSION(NPOINTS, 4) :: phis4                  ! Correction factor to B and vh
  
  REAL(SZ), DIMENSION(NPOINTS)    :: bs
  REAL(SZ), DIMENSION(NPOINTS, 4) :: bs4
  REAL(SZ), DIMENSION(NPOINTS)    :: vmBL
  REAL(SZ), DIMENSION(NPOINTS, 4) :: vmBL4
  INTEGER,  DIMENSION(NPOINTS, 4) :: quadFlag4
  REAL(SZ), DIMENSION(NPOINTS, 4) :: quadIR4
  REAL(SZ), DIMENSION(NQUADS)     :: vrQuadrant
  REAL(SZ), DIMENSION(NQUADS)     :: radius                 ! Wind radii - the distance

  INTEGER                         :: quad                   ! Quadrant counter

  REAL(SZ)                        :: latestRMax             ! most recently calculated value of fitted rmax
  REAL(SZ)                        :: latestAngle            ! angle of the most recently calculated node w.r.t. the storm location
  LOGICAL                         :: usequadrantVR
  LOGICAL                         :: useVMaxesBL


  CONTAINS


  !----------------------------------------------------------------
  ! S U B R O U T I N E   C A L C  I N T E N S I T Y  C H A N G E
  !----------------------------------------------------------------
  !> This subroutine calculates the intensity time change of a variable
  !> using second order mumerical accuracy and uneven spacing.
  !>
  !> On input:
  !>    var         The input variable (vector)
  !>    times       Time values (vector) at the center locations
  !>    order       The accuracy order required for the calculations (1, 2)
  !>                  <= 1: first order approximation for finite differences
  !>                  >= 2: second order approximation for finite differences
  !>
  !> On output:
  !>   calcInt      the calculated intensity change (df/dt)
  !>    status      error status (0 means no error)
  !----------------------------------------------------------------
  SUBROUTINE CalcIntensityChange(var, times, calcInt, status, order)

    USE PaHM_Global, ONLY : DEG2RAD
    USE PaHM_Utilities, ONLY : SphericalDistance

    IMPLICIT NONE

    REAL(SZ), DIMENSION(:), INTENT(IN)  :: var, times
    INTEGER, OPTIONAL, INTENT(IN)       :: order

    REAL(SZ), DIMENSION(:), INTENT(OUT) :: calcInt
    INTEGER, INTENT(OUT)                :: status

    INTEGER                             :: ordAcur
    REAL(SZ)                            :: dt1, dt2
    LOGICAL                             :: dt1OK, dt2OK
    REAL(SZ)                            :: val1, val2
    INTEGER                             :: iCnt, maxCnt

    status = 0
    maxCnt = 0

    CALL SetMessageSource("CalcIntensityChange")
    
    IF ((SIZE(SHAPE(var)) /= 1) .OR. (SIZE(SHAPE(times)) /= 1)) THEN
      WRITE(scratchMessage, '(a)') 'The rank of arrays var and times should be equal to 1 (vectors)'
      CALL AllMessage(ERROR, scratchMessage)

      CALL UnsetMessageSource()

      status = 1

      RETURN
    ELSE
      maxCnt = SIZE(var)
    END IF

    ordAcur = 2
    IF (PRESENT(order)) THEN
      IF (order <= 1) ordAcur = 1 
      IF (order  > 1) ordAcur = 2
    END IF
    IF (SIZE(var) < 3) ordAcur = 1

    ! Case 1st orded accuracy using backward differences
    IF (ordAcur == 1 )THEN
      DO iCnt = 2, maxCnt
        dt1 = times(iCnt) - times(iCnt - 1)
          dt1OK = (CompareReals(dt1, 0.0_SZ) /=0)

        val1 = 0.0_SZ
        IF (dt1OK) val1 = (var(iCnt) - var(iCnt - 1)) / dt1

        calcInt(iCnt) = val1
      END DO
      calcInt(1) = calcInt(2)

      CALL UnsetMessageSource()

      RETURN
    END IF

    ! Case 2nd order accuracy using Forward differences for the first point,
    ! backward differences for the last point and central differences in
    ! between points. Temporal spacing assumed to be uneven (general case).
    ! Forward, backward and central differences are all 2nd order accurate
    ! approximations.

    !----- Forward differences (first point)
    iCnt = 1
    dt1 = times(iCnt + 1) - times(iCnt)
      dt1OK = (CompareReals(dt1, 0.0_SZ) /=0)
    dt2 = times(iCnt + 2) - times(iCnt + 1)
      dt2OK = (CompareReals(dt2, 0.0_SZ) /=0)

    val1 = 0.0_SZ
    IF (dt1OK) val1 = (var(iCnt + 1) -  var(iCnt)) / dt1

    val2 = 0.0_SZ
    IF (dt2OK) val2 = (var(iCnt + 2) -  var(iCnt + 1)) / dt2

    IF (dt1OK .AND. dt2OK) THEN
      calcInt(iCnt) = ((2.0_SZ * dt1 + dt2) / (dt1 + dt2)) * val1 - (dt1 / (dt1 + dt2)) * val2
    ELSE IF (.NOT. dt1OK) THEN
      calcInt(iCnt) = val1
    ELSE
      calcInt(iCnt) = 2.0_SZ * val1 - val2
    END IF
    !----- Forward differences (first point)

    !----- Central differences
    DO iCnt = 2, maxCnt - 1
      ! Forward
      dt1 = times(iCnt + 1) - times(iCnt)
        dt1OK = (CompareReals(dt1, 0.0_SZ) /=0)
      ! Backward
      dt2 = times(iCnt) - times(iCnt - 1)
        dt2OK = (CompareReals(dt2, 0.0_SZ) /=0)

      val1 = 0.0_SZ
      IF (dt1OK) val1 = (var(iCnt + 1) -  var(iCnt)) / dt1

      val2 = 0.0_SZ
      IF (dt2OK) val2 = (var(iCnt) -  var(iCnt - 1)) / dt2

      IF (dt1OK .AND. dt2OK) THEN
        calcInt(iCnt) = (dt2 / (dt1 + dt2)) * val1 + (dt1 / (dt1 + dt2)) * val2
      ELSE IF (.NOT. dt1OK) THEN
        calcInt(iCnt) = val1
      ELSE
        calcInt(iCnt) = val2
      END IF
    END DO
    !----- Central differences

    !----- Backward differences (last point)
    iCnt = maxCnt
    dt1 = times(iCnt) - times(iCnt - 1)
      dt1OK = (CompareReals(dt1, 0.0_SZ) /=0)
    dt2 = times(iCnt - 1) - times(iCnt - 2)
      dt2OK = (CompareReals(dt2, 0.0_SZ) /=0)

    val1 = 0.0_SZ
    IF (dt1OK) val1 = (var(iCnt) -  var(iCnt - 1)) / dt1

    val2 = 0.0_SZ
    IF (dt2OK) val2 = (var(iCnt - 1) -  var(iCnt - 2)) / dt2

    IF (dt1OK .AND. dt2OK) THEN
      calcInt(iCnt) = ((2.0_SZ * dt1 + dt2) / (dt1 + dt2)) * val1 - (dt1 / (dt1 + dt2)) * val2
    ELSE IF (.NOT. dt1OK) THEN
      calcInt(iCnt) = val1
    ELSE
      calcInt(iCnt) = 2.0_SZ * val1 - val2
    END IF
    !----- Backward differences (last point)

    CALL UnsetMessageSource()

  END SUBROUTINE CalcIntensityChange

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   U V  T R A N S
  !----------------------------------------------------------------
  !> This subroutine calculates the translational velocity of a moving hurricane
  !> using second order mumerical accuracy and uneven spacing.
  !>
  !> On input:
  !>    lat         Latitude values (vector) of the center (degrees north)
  !>    lon         Longitude values (vector) of the center (degrees east )
  !>    times       Time values (vector) at the center locations (seconds)
  !>    order       The accuracy order required for the calculations (1, 2)
  !>                  <= 1: first order approximation for finite differences
  !>                  >= 2: second order approximation for finite differences
  !>
  !> On output:
  !>    u           x component of the translational velocities (m/s)
  !>    v           y component of the translational velocities (m/s)
  !>    status      error status (0 means no error)
  !----------------------------------------------------------------
  SUBROUTINE UVTrans(lat, lon, times, u, v, status, order)

    USE PaHM_Global, ONLY : DEG2RAD
    USE PaHM_Utilities, ONLY : SphericalDistance

    IMPLICIT NONE

    REAL(SZ), DIMENSION(:), INTENT(IN)  :: lat, lon, times
    INTEGER, OPTIONAL, INTENT(IN)       :: order

    REAL(SZ), DIMENSION(:), INTENT(OUT) :: u, v
    INTEGER, INTENT(OUT)                :: status

    INTEGER                             :: ordAcur
    REAL(SZ)                            :: dx1, dy1, dx2, dy2
    REAL(SZ)                            :: dt1, dt2
    LOGICAL                             :: dt1OK, dt2OK
    REAL(SZ)                            :: u1, u2, v1, v2
    INTEGER                             :: iCnt, maxCnt

    status = 0
    maxCnt = 0

    CALL SetMessageSource("UVTrans")

    IF ((SIZE(SHAPE(lat)) /= 1) .OR. (SIZE(SHAPE(lon)) /= 1) .OR. (SIZE(SHAPE(times)) /= 1)) THEN
      WRITE(scratchMessage, '(a)') 'The rank of arrays lat, lon and times should be equal to 1 (vectors)'
      CALL AllMessage(ERROR, scratchMessage)

      CALL UnsetMessageSource()

      status = 1

      RETURN
    ELSE
      maxCnt = SIZE(lat)
    END IF

    ordAcur = 2
    IF (PRESENT(order)) THEN
      IF (order <= 1) ordAcur = 1 
      IF (order  > 1) ordAcur = 2
    END IF
    IF (SIZE(lat) < 3) ordAcur = 1

    ! Case 1st orded accuracy using backward differences
    IF (ordAcur == 1 )THEN
      DO iCnt = 2, maxCnt
        dx1 = SphericalDistance(lat(iCnt - 1), lon(iCnt - 1), lat(iCnt - 1), lon(iCnt))
        dy1 = SphericalDistance(lat(iCnt - 1), lon(iCnt - 1), lat(iCnt), lon(iCnt - 1))
        dt1 = ABS(times(iCnt) - times(iCnt - 1))
          dt1OK = (CompareReals(dt1, 0.0_SZ) /=0)

        u1 = 0.0_SZ
        v1 = 0.0_SZ
        IF (dt1OK) THEN
          u1 = SIGN(dx1 / dt1, (lon(iCnt) - lon(iCnt - 1)))
          v1 = SIGN(dy1 / dt1, (lat(iCnt) - lat(iCnt - 1)))
        END IF

        u(iCnt) = u1
        v(iCnt) = v1
      END DO
      u(1) = u(2)
      v(1) = v(2)

      CALL UnsetMessageSource()

      RETURN
    END IF

    ! Case 2nd order accuracy using Forward differences for the first point,
    ! backward differences for the last point and central differences in
    ! between points. Temporal spacing assumed to be uneven (general case).
    ! Forward, backward and central differences are all 2nd order accurate
    ! approximations.

    !----- Forward differences (first point)
    iCnt = 1
    dx1 = SphericalDistance(lat(iCnt), lon(iCnt), lat(iCnt), lon(iCnt + 1))
    dy1 = SphericalDistance(lat(iCnt), lon(iCnt), lat(iCnt + 1), lon(iCnt))
    dt1 = ABS(times(iCnt + 1) - times(iCnt))
      dt1OK = (CompareReals(dt1, 0.0_SZ) /=0)
          
    dx2 = SphericalDistance(lat(iCnt + 1), lon(iCnt + 1), lat(iCnt + 1), lon(iCnt + 2))
    dy2 = SphericalDistance(lat(iCnt + 1), lon(iCnt + 1), lat(iCnt + 2), lon(iCnt + 1))
    dt2 = ABS(times(iCnt + 2) - times(iCnt + 1))
      dt2OK = (CompareReals(dt2, 0.0_SZ) /=0)

    u1 = 0.0_SZ
    v1 = 0.0_SZ
    IF (dt1OK) THEN
      u1 = SIGN(dx1 / dt1, (lon(iCnt + 1) - lon(iCnt)))
      v1 = SIGN(dy1 / dt1, (lat(iCnt + 1) - lat(iCnt)))
    END IF

    u2 = 0.0_SZ
    v2 = 0.0_SZ
    IF (dt2OK) THEN
      u2 = SIGN(dx2 / dt2, (lon(iCnt + 2) - lon(iCnt + 1)))
      v2 = SIGN(dy2 / dt2, (lat(iCnt + 2) - lat(iCnt + 1)))
    END IF

    IF (dt1OK .AND. dt2OK) THEN
      u(iCnt) = ((2.0_SZ * dt1 + dt2) / (dt1 + dt2)) * u1 - (dt1 / (dt1 + dt2)) * u2
      v(iCnt) = ((2.0_SZ * dt1 + dt2) / (dt1 + dt2)) * v1 - (dt1 / (dt1 + dt2)) * v2
    ELSE IF (.NOT. dt1OK) THEN
      u(iCnt) = u1
      v(iCnt) = v1
    ELSE
      u(iCnt) = 2.0_SZ * u1 - u2
      v(iCnt) = 2.0_SZ * v1 - v2
    END IF
    !----- Forward differences (first point)

    !----- Central differences
    DO iCnt = 2, maxCnt - 1
      ! Forward
      dx1 = SphericalDistance(lat(iCnt), lon(iCnt), lat(iCnt), lon(iCnt + 1))
      dy1 = SphericalDistance(lat(iCnt), lon(iCnt), lat(iCnt + 1), lon(iCnt))
      dt1 = ABS(times(iCnt + 1) - times(iCnt))
        dt1OK = (CompareReals(dt1, 0.0_SZ) /=0)
      ! Backward
      dx2 = SphericalDistance(lat(iCnt - 1), lon(iCnt - 1), lat(iCnt - 1), lon(iCnt))
      dy2 = SphericalDistance(lat(iCnt - 1), lon(iCnt - 1), lat(iCnt), lon(iCnt - 1))
      dt2 = ABS(times(iCnt) - times(iCnt - 1))
        dt2OK = (CompareReals(dt2, 0.0_SZ) /=0)

      u1 = 0.0_SZ
      v1 = 0.0_SZ
      IF (dt1OK) THEN
        u1 = SIGN(dx1 / dt1, (lon(iCnt + 1) - lon(iCnt)))
        v1 = SIGN(dy1 / dt1, (lat(iCnt + 1) - lat(iCnt)))
      END IF

      u2 = 0.0_SZ
      v2 = 0.0_SZ
      IF (dt2OK) THEN
        u2 = SIGN(dx2 / dt2, (lon(iCnt) - lon(iCnt - 1)))
        v2 = SIGN(dy2 / dt2, (lat(iCnt) - lat(iCnt - 1)))
      END IF

      IF (dt1OK .AND. dt2OK) THEN
        u(iCnt) = (dt2 / (dt1 + dt2)) * u1 + (dt1 / (dt1 + dt2)) * u2
        v(iCnt) = (dt2 / (dt1 + dt2)) * v1 + (dt1 / (dt1 + dt2)) * v2
      ELSE IF (.NOT. dt1OK) THEN
        u(iCnt) = u1
        v(iCnt) = v1
      ELSE
        u(iCnt) = u2
        v(iCnt) = v2
      END IF
    END DO
    !----- Central differences

    !----- Backward differences (last point)
    iCnt = maxCnt
    dx1 = SphericalDistance(lat(iCnt), lon(iCnt), lat(iCnt), lon(iCnt - 1))
    dy1 = SphericalDistance(lat(iCnt), lon(iCnt), lat(iCnt - 1), lon(iCnt))
    dt1 = ABS(times(iCnt) - times(iCnt - 1))
      dt1OK = (CompareReals(dt1, 0.0_SZ) /=0)

    dx2 = SphericalDistance(lat(iCnt - 1), lon(iCnt - 1), lat(iCnt - 1), lon(iCnt - 2))
    dy2 = SphericalDistance(lat(iCnt - 1), lon(iCnt - 1), lat(iCnt - 2), lon(iCnt - 1))
    dt2 = ABS(times(iCnt - 1) - times(iCnt - 2))
      dt2OK = (CompareReals(dt2, 0.0_SZ) /=0)

    u1 = 0.0_SZ
    v1 = 0.0_SZ
    IF (dt1OK) THEN
      u1 = SIGN(dx1 / dt1, (lon(iCnt) - lon(iCnt - 1)))
      v1 = SIGN(dy1 / dt1, (lat(iCnt) - lat(iCnt - 1)))
    END IF

    u2 = 0.0_SZ
    v2 = 0.0_SZ
    IF (dt2OK) THEN
      u2 = SIGN(dx2 / dt2, (lon(iCnt - 1) - lon(iCnt - 2)))
      v2 = SIGN(dy2 / dt2, (lat(iCnt - 1) - lat(iCnt - 2)))
    END IF
    
    IF (dt1OK .AND. dt2OK) THEN
      u(iCnt) = ((2.0_SZ * dt1 + dt2) / (dt1 + dt2)) * u1 - (dt1 / (dt1 + dt2)) * u2
      v(iCnt) = ((2.0_SZ * dt1 + dt2) / (dt1 + dt2)) * v1 - (dt1 / (dt1 + dt2)) * v2
    ELSE IF (.NOT. dt1OK) THEN
      u(iCnt) = u1
      v(iCnt) = v1
    ELSE
      u(iCnt) = 2.0_SZ * u1 - u2
      v(iCnt) = 2.0_SZ * v1 - v2
    END IF
    !----- Backward differences (last point)

    CALL UnsetMessageSource()

  END SUBROUTINE UVTrans

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   U V  T R A N S  P O I N T
  !----------------------------------------------------------------
  !> This subroutine calculates the translational velocity of a moving hurricane.
  !>
  !> On input:
  !>    lat1      Previous latitude  of center (degrees north)
  !>    lon1      Previous longitude of center (degrees east )
  !>    lat2      Current  latitude  of center (degrees north)
  !>    lon2      Current  longitude of center (degrees east)
  !>    time1     Previous time (seconds)
  !>    time1     Current  time (seconds)
  !>
  !> On output:
  !>    u         x component of translational velocity (m/s)
  !>    v         y component of translational velocity (m/s)
  !----------------------------------------------------------------
  SUBROUTINE UVTransPoint(lat1, lon1, lat2, lon2, time1, time2, u, v)

    USE PaHM_Global, ONLY : DEG2RAD
    USE PaHM_Utilities, ONLY : SphericalDistance

    IMPLICIT NONE

    ! Global variables
    REAL(SZ), INTENT(IN)  :: lat1, lon1, lat2, lon2
    REAL(SZ), INTENT(IN)  :: time1, time2
    REAL(SZ), INTENT(OUT) :: u, v

    ! Local variables
    REAL(SZ) :: dx, dy, dt
    LOGICAL  :: dtOK

    dx = SphericalDistance(lat1, lon1, lat1, lon2)
    dy = SphericalDistance(lat1, lon1, lat2, lon1)
    dt = ABS(time2 - time1)
      dtOK = (CompareReals(dt, 0.0_SZ) /=0)

    u = 0.0_SZ
    v = 0.0_SZ
    IF (dtOK) THEN
      u = SIGN(dx / dt, (lon2 - lon1))
      v = SIGN(dy / dt, (lat2 - lat1))
    END IF

  END SUBROUTINE UVTransPoint

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   N E W  V O R T E X
  !----------------------------------------------------------------
  !
  !> Create a new Vortex object.
  !>
  !> On input:
  !>    pn           Ambient surface pressure (mb)
  !>    pc           Surface pressure at center of storm (mb)
  !>    cLat         Latitude  of storm center (degrees north)
  !>    cLon         Longitude of storm center (degrees east )
  !>    vMax         Max sustained wind velocity in storm (knots)
  !>
  !> On output:
  !>    A new vortex is created with essential parameters calculated.
  !----------------------------------------------------------------
  SUBROUTINE NewVortex(pinf, p0, lat, lon, vm)

    USE PaHM_Global, ONLY : rhoAir, DEG2RAD, OMEGA, MB2PA, KT2MS

    IMPLICIT NONE

    REAL(SZ), INTENT(IN) :: pinf
    REAL(SZ), INTENT(IN) :: p0
    REAL(SZ), INTENT(IN) :: lat
    REAL(SZ), INTENT(IN) :: lon
    REAL(SZ), INTENT(IN) :: vm

    ! set instance variables
    pn = pinf
    pc = p0
    cLat = lat
    cLon = lon
    vMax = vm
!PV Check conversions
    ! evaluate basic physical params
    corio = 2.0_SZ * OMEGA * SIN(DEG2RAD * cLat)
    B = (vMax * KT2MS)**2 * rhoAir * EXP(1.0_SZ) / ((pn - pc) * MB2PA)
    B = MAX(MIN(B, 2.0_SZ), 1.0_SZ) ! limit B to range 1.0->2.5
!PV Data already have been converted
    ! added for compatibility of CalcRMaxes to use with simplified nws20
    bs(1:6) = B
    vmBL(1:6) = vMax

  END SUBROUTINE NewVortex

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   N E W  V O R T E X  F U L L
  !----------------------------------------------------------------
  !> A new vortex is created for the full gradient wind balance
  !>
  !----------------------------------------------------------------
  SUBROUTINE NewVortexFull(pinf, p0, lat, lon, vm)

    USE PaHM_Global, ONLY : rhoAir, DEG2RAD, KT2MS, OMEGA, MB2PA

    IMPLICIT NONE

    REAL(SZ), INTENT(IN) :: pinf
    REAL(SZ), INTENT(IN) :: p0
    REAL(SZ), INTENT(IN) :: lat
    REAL(SZ), INTENT(IN) :: lon
    REAL(SZ), INTENT(IN) :: vm
 
    ! set instance variables
    pn   = pinf
    pc   = p0
    cLat = lat
    cLon = lon
    vMax = vm
 
    ! evaluate basic physical params
    corio = 2.0_SZ * OMEGA * SIN(DEG2RAD * cLat)
    B = (vMax * KT2MS)**2 * rhoAir * EXP(1.0_SZ) / ((pn - pc) * MB2PA)
    phi       = 1.0_SZ
    bs(1:6)   = B
    phis(1:6) = phi
    vmBL(1:6) = vMax

    ! Jie 2013.01 
    ! B = MAX(MIN(B, 2.0_SZ), 1.0_SZ) ! limit B to range 1.0->2.5

  END SUBROUTINE NewVortexFull

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   S E T  V O R T E X
  !----------------------------------------------------------------
  !> Set basic parameter for a new Vortex object.
  !>
  !> On input:
  !>    pinf         hurricane Ambient pressure
  !>    p0           hurricane central pressure
  !>    lat          Latitude  of storm center (degrees north)
  !>    lon          Longitude of storm center (degrees east )
  !>
  !> On output:
  !>    Aim is to define pn, pc, and corio
  !----------------------------------------------------------------
  SUBROUTINE SetVortex(pinf, p0, lat, lon)

    USE PaHM_Global, ONLY : DEG2RAD, OMEGA

    IMPLICIT NONE

    REAL(SZ), INTENT(IN) :: pinf
    REAL(SZ), INTENT(IN) :: p0
    REAL(SZ), INTENT(IN) :: lat
    REAL(SZ), INTENT(IN) :: lon

    ! set instance variables
    pn   = pinf
    pc   = p0
    cLat = lat
    cLon = lon

    ! evaluate basic physical params
    corio = 2.0_SZ * OMEGA * SIN(DEG2RAD * cLat)

  END SUBROUTINE SetVortex

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   S E T   R M A X E S
  !----------------------------------------------------------------
  SUBROUTINE SetRMaxes(rMaxW)

    IMPLICIT NONE

    REAL(SZ), DIMENSION(4), INTENT(IN) :: rMaxW
    INTEGER :: i

    DO i = 1, 4
       rMaxes(i + 1) = rMaxW(i)
    END DO

  END SUBROUTINE SetRMaxes

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   G E T  R M A X E S
  !----------------------------------------------------------------
  SUBROUTINE GetRMaxes(rMaxW)

    IMPLICIT NONE

    REAL(SZ), DIMENSION(4), INTENT(OUT) :: rMaxW

    INTEGER :: i

    DO i = 1, 4
      rMaxW(i) = rMaxes(i + 1)
    END DO

  END SUBROUTINE GetRMaxes

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   C A L C  R M A X E S
  !----------------------------------------------------------------
  !> Calculate the radius of maximum winds for all storm quadrants.
  !>
  !> On input:
  !>    none
  !>
  !> On output:
  !>    rMax    radius of maximum winds (nm) in all quadrants, plus
  !>            2 extra values to tie down circular periodicity
  !  Jie 2014.07 Modified with quadrant-varying vmBL, which not only
  !  works for nws19 but for the simplified nws20
  !----------------------------------------------------------------
  SUBROUTINE CalcRMaxes()

    IMPLICIT NONE

    REAL(SZ)            :: root        ! Radius of maximum winds
    REAL(SZ), PARAMETER :: INNERRADIUS = 1.0_SZ
    REAL(SZ), PARAMETER :: OUTERRADIUS = 400.0_SZ
    REAL(SZ), PARAMETER :: ACCURACY    = 0.0001_SZ
    REAL(SZ), PARAMETER :: ZOOM        = 0.01_SZ
    INTEGER , PARAMETER :: ITERMAX     = 3
    REAL(SZ)            :: r1, r2, r3, r4, dr
    REAL(SZ)            :: vicinity
    INTEGER             :: n, iter

    !-----------------------------
    ! Loop over quadrants of storm
    !-----------------------------
    DO n = 1, nQuads
      ! set B and vMax values for each quadrant
      ! for nws19, B and vMax are constant
      ! for simplified nws20, B is constant, while vMax is not
      B = bs(n + 1)
      vMax = vmBL(n + 1)

      quad = n
      root = -1.0_SZ
      r1   = INNERRADIUS
      r2   = OUTERRADIUS
      dr   = 1.0_SZ
      DO iter = 1, ITERMAX
        root = FindRoot(VhWithCori, r1, r2, dr, r3, r4)
        r1   = r3
        r2   = r4
        dr   = dr * ZOOM
      END DO

      ! determine if rMax is actually in the vicinity of the
      ! isotach radius that we are using to solve for rMax,
      ! and if so, take another shot at finding the
      ! rMax using the gradient wind balance that neglects
      ! coriolis (and is appropriate in the vicinity of rMax)
      vicinity = ABS(root - radius(quad)) / root
      IF ((root < 0.0_SZ) .OR. (vicinity <= 0.0_SZ)) THEN
        r1 = INNERRADIUS
        r2 = OUTERRADIUS
        dr = 1.0_SZ
        DO iter = 1, ITERMAX
          root = FindRoot(VhNoCori, r1, r2, dr, r3, r4)
          r1 = r3
          r2 = r4
          dr = dr * ZOOM
        END DO
      END IF

      rMaxes(n + 1) = root
    END DO

  END SUBROUTINE CalcRMaxes

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   C A L C  R M A X E S  F U L L
  !----------------------------------------------------------------
  !> Calculate the radius of maximum winds for all storm quadrants.
  !> Solving the full gradient wind equation without the assumption
  !> of cyclostrohpic balance.
  !>
  !> On input:
  !>    none
  !>
  !> On output:
  !>    rMax    radius of maximum winds (nm) in all quadrants, plus
  !>            2 extra values to tie down circular periodicity
  !>
  ! Jie 2013.02 added looping procedures to calculate bs and phis
  !----------------------------------------------------------------
  SUBROUTINE CalcRMaxesFull()

    USE PaHM_Global, ONLY : rhoAir, NM2M, KT2MS, MB2PA

    IMPLICIT NONE

    REAL(SZ)            :: root           ! Radius of maximum winds
    REAL(SZ), PARAMETER :: INNERRADIUS = 1.0_SZ
    REAL(SZ), PARAMETER :: OUTERRADIUS = 500.0_SZ
    REAL(SZ), PARAMETER :: ACCURACY    = 0.0001_SZ
    REAL(SZ), PARAMETER :: ZOOM        = 0.01_SZ
    INTEGER , PARAMETER :: ITERMAX     = 3
    REAL(SZ)            :: r1, r2, r3, r4, dr
    INTEGER             :: n, iter, noRootFlag
    REAL(SZ)            :: bNew, bNew1
    REAL(SZ)            :: phiNew
    INTEGER, PARAMETER  :: cont = 400     ! Max # of iterations
    INTEGER             :: iCont, ibCont  ! iteration counter

    211 FORMAT(a7, x ,i2, x, a38)

    !-----------------------------
    ! Loop over quadrants of storm
    !-----------------------------
    DO n = 1, nQuads
      noRootFlag = 0

      ! initialize B and phi values for each quadrant
      B    = bs(n + 1)
      phi  = phis(n + 1)
      vMax = vmBL(n + 1)

      ! Loop the root-solving process to converge B, for in the
      ! new wind formulation, B is a function of rMax, vMax, f, and phi 
      DO iCont = 1, cont ! logical expre. is at the end to exit the loop
        noRootFlag = 0
        quad = n
        root = -1.0_SZ
        r1 = INNERRADIUS
        r2 = OUTERRADIUS
        dr = 1.0_SZ
        DO iter = 1, ITERMAX
          root = FindRoot(VhWithCoriFull, r1, r2, dr, r3, r4)
          r1 = r3
          r2 = r4
          dr = dr * ZOOM
        END DO

        ! avoid invalid B value when root is not found        
        IF (root < 0.0_SZ) THEN
        !  r1 = INNERRADIUS
        !  r2 = OUTERRADIUS
        !  dr = 1.0_SZ
        !  DO iter = 1, ITERMAX
        !    root = FindRoot(VhNoCori, r1, r2, dr, r3, r4)
        !    r1 = r3
        !    r2 = r4
        !    dr = dr * ZOOM
        !  END DO
          root = 1.0 * radius(quad)
          noRootFlag = 1
        END IF

        rMaxes(n + 1) = root  
              
        ! Jie 2013.02
        ! determine if B converges, if yes, break loop and assign 
        ! values to rMaxes, if not, continue the loop to re-calculate 
        ! root and re-evaluate bs
        phiNew = 1 + vMax * KT2MS * root * NM2M * corio /                                &
                (B * ((vMax * KT2MS)**2 + vMax * KT2MS * root * NM2M * corio))
        bNew   = ((vMax * KT2MS)**2 + vMax * KT2MS * root * NM2M * corio) *              &
                 rhoAir * EXP(phiNew) / (phiNew * (pn - pc) * MB2PA)
        DO ibCont = 1, cont
          bNew1 = bNew
          phiNew = 1 + vMax * KT2MS * root * NM2M * corio /                              &
                   (bNew * ((vMax * KT2MS)**2 + vMax * KT2MS * root * NM2M * corio))
          bNew   = ((vMax * KT2MS)**2 + vMax * KT2MS * root * NM2M * corio) *            &
                   rhoAir * EXP(phiNew) / (phiNew * (pn - pc) * MB2PA)

          IF (ABS(bNew - bNew1) <= 0.01_SZ) EXIT
        END DO

        ! debug with aswip          
        !IF (ibCont >= cont) THEN
        !  WRITE(1111, 211) "iquad=", n, "bNew did not fully converge, procede"
        !END IF
        ! end debug with aswip         

        IF (ABS(B - bNew) <= 0.01_SZ) EXIT
 
        ! update B and phi for next iteration
        ! warning: modifications made here also affect other subroutines
        B   = bNew
        phi = phiNew
      END DO !iCont = 1, cont 

      ! update to the latest values for aswip output
      bs(n + 1)   = bNew
      phis(n +1 ) = phiNew

      ! debug with aswip         
      !IF (iCont >= cont) THEN
      !  WRITE(1111, 211) "iquad=", n, "B did not fully converge, procede"
      !END IF
      ! end debug with aswip 

      ! determine if rMax is actually in the vicinity of the
      ! isotach radius that we are using to solve for rMax,
      ! and if so, take another shot at finding the
      ! rMax using the gradient wind equation that neglects
      ! coriolis (and is appropriate in the vicinity of rMax)
      ! Jie 2013.01
      !vicinity = ABS(root - radius(quad)) / root
      IF (noRootFlag == 1) THEN
         WRITE(*, *) "iquad=", n, "No root found, return dist. to Isotach"  
      END IF
    END DO  !n = 1, nQuads

  END SUBROUTINE CalcRMaxesFull

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   F I T   R M A X E S
  !----------------------------------------------------------------
  !  RJW 07 - 2009
  !> Calculates the coefficients that fit the given
  !> radius of maximum winds for all storm quadrants.
  !>
  !> On input:
  !>   rMax in all 4 quadrants plus 2 extra values to tie down circular periodicity
  !>
  !> On output:
  !>    rMax    radius of maximum winds (nm) in all quadrants, plus
  !>            2 extra values to tie down circular periodicity
  !----------------------------------------------------------------
  SUBROUTINE FitRMaxes()

    IMPLICIT NONE

    ! Generate 2 additional (theta, rMax) points for curve-fit
    rMaxes(1) = rMaxes(5)
    rMaxes(6) = rMaxes(2)

  END SUBROUTINE FitRMaxes

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   F I T   R M A X E S  4
  !----------------------------------------------------------------
  SUBROUTINE FitRMaxes4()

    IMPLICIT NONE

    ! Generate 2 additional points for curve-fit
    quadFlag4(1, 1:4) = quadFlag4(5, 1:4)
    quadFlag4(6, 1:4) = quadFlag4(2, 1:4)
    
    quadIR4(1, 1:4) = quadIR4(5, 1:4)
    quadIR4(6, 1:4) = quadIR4(2, 1:4)
          
    rMaxes4(1, 1:4) = rMaxes4(5, 1:4)
    rMaxes4(6, 1:4) = rMaxes4(2, 1:4)
    
    bs4(1, 1:4) = bs4(5, 1:4)
    bs4(6, 1:4) = bs4(2, 1:4)
    
    phis4(1, 1:4) = phis4(5, 1:4)
    phis4(6, 1:4) = phis4(2, 1:4)  
    
    vmBL4(1, 1:4) = vmBL4(5, 1:4)
    vmBL4(6, 1:4) = vmBL4(2, 1:4) 
 
  END SUBROUTINE FitRMaxes4

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   S E T  V M A X E S  B L
  !----------------------------------------------------------------
  SUBROUTINE SetVMaxesBL(vMaxW)

    IMPLICIT NONE

    REAL(SZ), DIMENSION(4), INTENT(IN) :: vMaxW

    INTEGER :: i

    DO i = 1, 4
      vmBL(i + 1) = vMaxW(i)
    END DO

  END SUBROUTINE SetVMaxesBL

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   G E T  V M A X E S  B L
  !----------------------------------------------------------------
  SUBROUTINE GetVMaxesBL(vMaxW)

    IMPLICIT NONE

    REAL(SZ), DIMENSION(4), INTENT(OUT) :: vMaxW

    INTEGER :: i

    DO i = 1, 4
      vMaxW(i) = vmBL(i + 1)
    END DO

  END SUBROUTINE GetVMaxesBL

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   S E T  U S E  V M A X E S  B L
  !----------------------------------------------------------------
  SUBROUTINE SetUseVMaxesBL(u)

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: u

    useVMaxesBL = u

  END SUBROUTINE SetUseVMaxesBL

  !----------------------------------------------------------------
  ! S U B R O U T I N E   S E T  S H A P E  P A R A M E T E R
  !----------------------------------------------------------------
  SUBROUTINE SetShapeParameter(param)

    IMPLICIT NONE

    REAL(SZ) :: param

    B = param

  END SUBROUTINE SetShapeParameter

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   G E T  S H A P E  P A R A M E T E R
  !----------------------------------------------------------------
  REAL(SZ) FUNCTION GetShapeParameter() RESULT(myValOut)

    IMPLICIT NONE

    myValOut = B

    RETURN

  END FUNCTION GetShapeParameter

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   G E T  S H A P E  P A R A M E T E R S
  !----------------------------------------------------------------
  FUNCTION GetShapeParameters() RESULT(myValOut)

    IMPLICIT NONE

    REAL(SZ), DIMENSION(4) :: myValOut

    INTEGER :: i

    DO i = 1, 4
      myValOut(i) = bs(i + 1)
    END DO

    RETURN

  END FUNCTION GetShapeParameters

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   G E T  P H I  F A C T O R S
  !----------------------------------------------------------------
  FUNCTION GetPhiFactors() RESULT(myValOut)

    IMPLICIT NONE

    REAL(SZ), DIMENSION(4) :: myValOut

    INTEGER :: i

    DO i = 1, 4
      myValOut(i) = phis(i + 1)
    END DO

    RETURN

  END FUNCTION GetPhiFactors

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   S E T  I S O T A C H  R A D I I
  !----------------------------------------------------------------
  SUBROUTINE SetIsotachRadii(ir)

    IMPLICIT NONE

    REAL(SZ), DIMENSION(4), INTENT(IN) :: ir

    radius(:) = ir(:)

  END SUBROUTINE SetIsotachRadii

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   S E T  I S O T A C H  W I N D  S P E E D S
  !----------------------------------------------------------------
  SUBROUTINE SetIsotachWindSpeeds(vrQ)

    IMPLICIT NONE

    REAL(SZ), DIMENSION(4), INTENT(IN) :: vrQ

    vrQuadrant(:) = vrQ(:)

  END SUBROUTINE SetIsotachWindSpeeds

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   S E T  U S E  Q U A D R A N T  V R
  !----------------------------------------------------------------
  SUBROUTINE SetUsequadrantVR(u)

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: u

    usequadrantVR = u

  END SUBROUTINE SetUsequadrantVR

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   G E T  L A T E S T  A N G L E
  !----------------------------------------------------------------
  LOGICAL FUNCTION GetUsequadrantVR() RESULT(myValOut)

    IMPLICIT NONE

    myValOut = usequadrantVR

  END FUNCTION GetUsequadrantVR

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   S P I N T E R P
  !----------------------------------------------------------------
  !
  !> Spatial Interpolation function based on angle and r.
  !>
  !> On input:
  !>    angle        Azimuthal angle (degrees)
  !>    r            Distnace to storm Center (nm)
  !>
  !> On output:
  !>    interpolated value for rMax/vMax/B          
  !>
  !> INTEGER validIsot is used as a marker to indicate how many isotachs
  !> are available in a certain quadrant
  !> SELECT CASE(validIsot)
  !> CASE(1): 1 situation
  !> CASE(2): 3 situations
  !> CASE(3): 4 situations
  !> CASE(4): 5 situations                    
  !----------------------------------------------------------------
  REAL(SZ) FUNCTION SpInterp(angle, dist, opt) RESULT(myValOut)

    IMPLICIT NONE

    REAL(SZ), INTENT(IN)            :: angle, dist  
    INTEGER, INTENT(IN)             :: opt
    REAL(SZ), DIMENSION(NPOINTS, 4) :: param 
    REAL(SZ)                        :: temp1, temp2  
    REAL(SZ)                        :: deltaAngle
    INTEGER                         :: iQuad

    IF (opt == 1) THEN
      param = rMaxes4
    ELSE IF (opt == 2) THEN
      param = bs4
    ELSE IF (opt == 3) THEN
      param = vmBL4    
    END IF   
    
    IF (angle <= 45.0_SZ) THEN
      iQuad = 5
      deltaAngle = 45.0_SZ + angle
    ELSE IF (angle <= 135.0_SZ) THEN
      iQuad = 2
      deltaAngle = angle - 45.0_SZ
    ELSE IF (angle <= 225.0_SZ) THEN
      iQuad = 3
      deltaAngle = angle - 135.0_SZ
    ELSE IF (angle <= 315.0_SZ) THEN
      iQuad = 4
      deltaAngle = angle - 225.0_SZ
    ELSE IF (angle > 315.0_SZ) THEN
      iQuad = 5
      deltaAngle = angle - 315.0_SZ
    END IF

    ! nearest neighbor weighted interpolation     
    IF ( deltaAngle < 1.0_SZ ) THEN
      myValOut = InterpR(param, iQuad, dist)
    ELSE IF (deltaAngle > 89.0_SZ) THEN
      myValOut = InterpR(param, iQuad + 1, dist)
    ELSE
      temp1 = InterpR(param, iQuad, dist)
      temp2 = InterpR(param, iQuad + 1, dist)
      myValOut = (temp1 / deltaAngle**2 + temp2 / (90.0 - deltaAngle)**2) /   &
              (1.0_SZ / deltaAngle**2 + 1.0_SZ / (90.0_SZ - deltaAngle)**2)
    END IF

  END FUNCTION SpInterp

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   I N T E R P R
  !----------------------------------------------------------------
  REAL(SZ) FUNCTION InterpR(quadVal, quadSel, quadDis) RESULT(myValOut)

    IMPLICIT NONE

    REAL(SZ), DIMENSION(NPOINTS, 4), INTENT(IN) :: quadVal
    INTEGER, INTENT(IN)                         :: quadSel
    REAL(SZ), INTENT(IN)                        :: quadDis

    REAL(SZ)                                    :: fac
    INTEGER                                     :: totalISot
    
    totalISot = SUM(quadFlag4(quadSel, :))
    SELECT CASE(totalISot)
      CASE(1)
        myValOut = quadVal(quadSel, MAXLOC(quadFlag4(quadSel, :), 1)) 
      CASE(2)
        IF (quadDis > quadIR4(quadSel, 1)) THEN
          myValOut = quadVal(quadSel, 1)
        ELSE IF (quadDis > quadIR4(quadSel, 2)) THEN
          fac = (quadDis - quadIR4(quadSel, 2)) / (quadIR4(quadSel, 1) - quadIR4(quadSel, 2))
          myValOut = quadVal(quadSel, 1) * fac + quadVal(quadSel, 2) * (1 - fac)
        ELSE
           myValOut = quadVal(quadSel, 2)
        END IF
      CASE(3)
        IF (quadDis > quadIR4(quadSel, 1)) THEN
          myValOut = quadVal(quadSel, 1)
        ELSE IF (quadDis > quadIR4(quadSel, 2)) THEN
          fac = (quadDis - quadIR4(quadSel, 2)) / (quadIR4(quadSel, 1) - quadIR4(quadSel, 2))
          myValOut = quadVal(quadSel, 1) * fac + quadVal(quadSel, 2) * (1 - fac)
        ELSE IF (quadDis > quadIR4(quadSel, 3)) THEN
          fac = (quadDis - quadIR4(quadSel, 3)) / (quadIR4(quadSel, 2) - quadIR4(quadSel, 3))
          myValOut = quadVal(quadSel, 2) * fac + quadVal(quadSel, 3) * (1 - fac)
        ELSE
          myValOut = quadVal(quadSel, 3)
        END IF
      CASE(4)
        IF (quadDis > quadIR4(quadSel, 1)) THEN
          myValOut = quadVal(quadSel, 1)
        ELSE IF (quadDis > quadIR4(quadSel, 2)) THEN
          fac = (quadDis - quadIR4(quadSel, 2)) / (quadIR4(quadSel, 1) - quadIR4(quadSel, 2))
          myValOut = quadVal(quadSel, 1) * fac + quadVal(quadSel, 2) * (1 - fac)
        ELSE IF (quadDis > quadIR4(quadSel, 3)) THEN
          fac = (quadDis - quadIR4(quadSel, 3)) / (quadIR4(quadSel, 2) - quadIR4(quadSel, 3))
          myValOut = quadVal(quadSel, 2) * fac + quadVal(quadSel, 3) * (1 - fac)
        ELSE IF (quadDis > quadIR4(quadSel, 4)) THEN
          fac = (quadDis - quadIR4(quadSel, 4)) / (quadIR4(quadSel, 3) - quadIR4(quadSel, 4))
          myValOut = quadVal(quadSel, 3) * fac + quadVal(quadSel, 4) * (1 - fac)
        ELSE
          myValOut = quadVal(quadSel, 4)
        END IF         
      CASE default
        ! For whatever reason if our algorithm fails, add the following
        ! line to avoid run-time errors
        myValOut = quadVal(quadSel, MAXLOC(quadFlag4(quadSel, :), 1)) 
        WRITE(*, *) "ERROR: InterpR failed in nws20get." !PV remove it of modify it?
    END SELECT
         
  END FUNCTION InterpR

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   R M W
  !----------------------------------------------------------------
  !
  !> Calculate the radius of maximum winds.
  !>
  !> On input:
  !>    angle        Azimuthal angle (degrees)
  !>
  !> On output:
  !>    Rmw          Radius of maximum winds (meters) from curve fit
  !>                      I DO NOT BELIEVE IT IS IN METERS rjw
  !----------------------------------------------------------------
  REAL(SZ) FUNCTION Rmw(angle) RESULT(myValOut)

    IMPLICIT NONE

    REAL(SZ), INTENT(IN) :: angle
    INTEGER              :: baseQuadrant
    REAL(SZ)             :: deltaAngle

    IF (angle <= 45.0_SZ) THEN
      baseQuadrant = 5
      deltaAngle = 45.0_SZ + angle
    ELSE IF (angle <= 135.0_SZ) THEN
      baseQuadrant = 2
      deltaAngle = angle - 45.0_SZ
    ELSE IF (angle <= 225.0_SZ) THEN
      baseQuadrant = 3
      deltaAngle = angle - 135.0_SZ
    ELSE IF (angle <= 315.0_SZ) THEN
      baseQuadrant = 4
      deltaAngle = angle - 225.0_SZ
    ELSE IF (angle > 315.0_SZ) THEN
      baseQuadrant = 5
      deltaAngle = angle - 315.0_SZ
    END IF

    ! nearest neighbor weighted interpolation
    IF ( deltaAngle < 1.0_SZ ) THEN
      myValOut = rMaxes(baseQuadrant)   ! avoid div by zero
    ELSE IF ( deltaAngle > 89.0_SZ ) THEN
      myValOut = rMaxes(baseQuadrant + 1) ! avoid div by zero
    ELSE
      myValOut = (rMaxes(baseQuadrant) / deltaAngle**2 +                       &
                  rMaxes(baseQuadrant + 1) / (90.0 - deltaAngle)**2) /         &
                 (1.0_SZ / deltaAngle**2 + 1.0_SZ / (90.0_SZ - deltaAngle)**2)
    END IF

    ! linearly interpolate
    !myValOut = (deltaAngle / 90.0_SZ) *                               &
    !           (rMaxes(baseQuadrant + 1) - rMaxes(baseQuadrant)) +   &
    !           rMaxes(baseQuadrant)

  END FUNCTION Rmw

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   U V P
  !----------------------------------------------------------------
  !
  !> Calculate (u, v) wind components and surface pressure from an
  !> asymmetric hurricane wind model.
  !>
  !> On input:
  !>    lat         Latitude  of nodal point (degrees north)
  !>    lon         Longitude of nodal point (degrees east )
  !>    uTrans      x component of translational velocity (kts)
  !>    vTrans      y component of translational velocity (kts)
  !>
  !> On output:
  !>    u           x component of wind velocity at nodal point (m/s)
  !>    v           y component of wind velocity at nodal point (m/s)
  !>    p           Surface pressure at nodal point (Pa)
  !>
  !> Internal parameters:
  !>    dampRadii   How far out (# of rMax radii) to begin damping
  !>                out the translational velocity
  !>
  !> Note:
  !>    Subroutine directly accesses global class instance variables
  !----------------------------------------------------------------
  SUBROUTINE UVP(lat, lon, uTrans, vTrans, u, v, p)

    USE PaHM_Global, ONLY : windReduction, ONE2TEN, DEG2RAD, RAD2DEG, MB2PA, KT2MS, NM2M, M2NM, REARTH

    IMPLICIT NONE

    REAL(SZ), INTENT(IN)  :: lat
    REAL(SZ), INTENT(IN)  :: lon
    REAL(SZ), INTENT(IN)  :: uTrans
    REAL(SZ), INTENT(IN)  :: vTrans

    REAL(SZ), INTENT(OUT) :: u
    REAL(SZ), INTENT(OUT) :: v
    REAL(SZ), INTENT(OUT) :: p

    REAL(SZ)              :: transSpdX  !NWS8-style translation speed
    REAL(SZ)              :: transSpdY  !NWS8-style translation speed

    REAL(SZ)              :: dx
    REAL(SZ)              :: dy
    REAL(SZ)              :: dist
    REAL(SZ)              :: rmx
    REAL(SZ)              :: angle
    REAL(SZ)              :: speed
    REAL(SZ)              :: uf
    REAL(SZ)              :: vf
    REAL(SZ)              :: percentCoriolis
    REAL(SZ)              :: speedAtRMax
    REAL(SZ)              :: vMaxFactor

    !------------------------------------------------------
    ! Calculate distance and angle between eye of hurricane
    ! and input nodal point
    !------------------------------------------------------
    dx = DEG2RAD * REARTH * (lon - cLon) * COS(DEG2RAD * cLat)
    dy = DEG2RAD * REARTH * (lat - cLat)
    dist = SQRT(dx * dx + dy * dy)

    !----------------------------------------
    ! Handle special case at eye of hurricane
    ! in eye velocity is zero not translational velocity
    !----------------------------------------
    IF (dist < 1.0_SZ) THEN
      u = 0.0_SZ
      v = 0.0_SZ
      p = pc * MB2PA

      RETURN
    END IF

    dist = M2NM * dist

    angle = 360.0_SZ + RAD2DEG * ATAN2(dx, dy)
    IF (angle > 360.0_SZ) angle = angle - 360.0_SZ
 
    latestAngle = angle
    rmx = Rmw(angle)
    latestRMax = rmx

    !---------------------------------------------------
    ! Compute (u,v) wind velocity components from the
    ! asymmetric hurricane vortex.
    !
    ! Note: the vortex winds are valid at the top of the
    ! surface layer, so reduce the winds to the surface.
    ! Also convert the winds from max sustained 1-minute
    ! averages to 10-minute averages for the storm surge
    ! model.
    !---------------------------------------------------
    percentCoriolis = 1.0_SZ
    speed = SQRT((vMax * KT2MS)**2 * (rmx / dist)**B * EXP(1.0_SZ - (rmx / dist)**B) +   &
                 (NM2M * dist * percentCoriolis * corio / 2.0_SZ)**2)                    &
                - NM2M * dist * percentCoriolis * corio / 2.0_SZ

    ! calculate the wind speed (m/s) at rMax, using
    ! equation that includes full coriolis
    speedAtRMax = SQRT((vMax * KT2MS)**2 * EXP(0.0_SZ) +                      &
                       (NM2M * dist * percentCoriolis * corio / 2.0_SZ)**2)   &
                      - NM2M * dist * percentCoriolis * corio / 2.0_SZ

    ! calculate a factor to place the velocity profile so that
    ! it hits vMax
    vMaxFactor = vMax * KT2MS / speedAtRMax

    ! jgf20111007: Calculate NWS8-like translation speed
    transSpdX = (ABS(speed / speedAtRMax)) * uTrans * KT2MS
    transSpdY = (ABS(speed / speedAtRMax)) * vTrans * KT2MS

    speed = speed * vMaxFactor

    ! now reduce the wind speed to the surface
    speed = speed * windReduction

    u = -speed * COS(DEG2RAD * angle)
    v =  speed * SIN(DEG2RAD * angle)

    ! Alter wind direction by adding a frictional inflow angle
    CALL Rotate(u, v, FAng(dist, rmx), cLat, uf, vf)
    u = uf
    v = vf
    !
    ! jgf20111007: Add in the translation velocity
    u = u + transSpdX
    v = v + transSpdY
    !
    ! convert from 1 minute averaged winds to 10 minute averaged
    ! winds for use in ADCIRC
    u = u * ONE2TEN
    v = v * ONE2TEN

    ! Compute surface pressure from asymmetric hurricane vortex
    p = MB2PA * (pc + (pn - pc) * EXP(-(rmx / dist)**B))

    ! cut off the vortex field after 401nm
    ! TODO: 401nm should be replaced with something less
    ! arbitrary ... and find a better way to blend this
    !IF ( dist > 401.0_SZ ) THEN
    !  u = 0.0_SZ
    !  v = 0.0_SZ
    !  p = MB2PA * pn
    !END IF

  END SUBROUTINE UVP

!================================================================================


  !----------------------------------------------------------------
  ! S U B R O U T I N E   U V P R
  !----------------------------------------------------------------
  !
  !> Calculate (u, v) wind components and surface pressure from an
  !> asymmetric hurricane wind model.
  !>
  !> On input:
  !>    pinf         hurricane Ambient pressure
  !>    p0           hurricane central pressure
  !>    iDist        dist to hurricane center in nautical mile
  !>    iRmx         Rmw
  !>    iAngle       Azimuth Angle
  !>    iB           Holland B parameter
  !>    iVm          vortex maximum velocity at upper boundary
  !>    iPhi         vortex correction factor
  !>    uTrans       x component of translational velocity (kts)
  !>    vTrans       y component of translational velocity (kts)
  !>
  !> On output:
  !>    u            x component of wind velocity at nodal point (m/s)
  !>    v            y component of wind velocity at nodal point (m/s)
  !>    p            Surface pressure at nodal point (Pa)
  !>
  !> Internal parameters:
  !>    dampRadii   How far out (# of rMax radii) to begin damping
  !>                out the translational velocity
  !>
  !> Note:
  !>    Subroutine directly accesses global class instance variables
  !----------------------------------------------------------------
  SUBROUTINE UVPR(iDist, iAngle, iRmx, iRmxTrue, iB, iVm, iPhi, &
                  uTrans, vTrans, geof, u, v, p)

    USE PaHM_Global, ONLY : windReduction, ONE2TEN, DEG2RAD, MB2PA, KT2MS, NM2M

    IMPLICIT NONE

    REAL(SZ), INTENT(IN)  :: iDist
    REAL(SZ), INTENT(IN)  :: iAngle
    REAL(SZ), INTENT(IN)  :: iRmx
    REAL(SZ), INTENT(IN)  :: iRmxTrue
    REAL(SZ), INTENT(IN)  :: iB
    REAL(SZ), INTENT(IN)  :: iVm
    REAL(SZ), INTENT(IN)  :: iPhi
    REAL(SZ), INTENT(IN)  :: uTrans
    REAL(SZ), INTENT(IN)  :: vTrans
    INTEGER , INTENT(IN)  :: geof

    REAL(SZ), INTENT(OUT) :: u
    REAL(SZ), INTENT(OUT) :: v
    REAL(SZ), INTENT(OUT) :: p

    REAL(SZ)              :: transSpdX  !NWS8-style translation speed
    REAL(SZ)              :: transSpdY  !NWS8-style translation speed
    REAL(SZ)              :: rmx
    REAL(SZ)              :: speed
    REAL(SZ)              :: uf
    REAL(SZ)              :: vf
    REAL(SZ)              :: percentCoriolis
      
    rmx  = iRmx
    B    = iB
    vMax = iVm
    phi  = iPhi

    !----------------------------------------
    ! Handle special case at eye of hurricane
    ! in eye velocity is zero not translational velocity
    !----------------------------------------
    IF (iDist < 1.0_SZ) THEN
      u = 0.0_SZ
      v = 0.0_SZ
      p = pc * MB2PA

      RETURN
    END IF

    !---------------------------------------------------
    ! Compute (u, v) wind velocity components from the
    ! asymmetric hurricane vortex.
    !
    ! Note: the vortex winds are valid at the top of the
    ! surface layer, so reduce the winds to the surface.
    ! Also convert the winds from max sustained 1-minute
    ! averages to 10-minute averages for the storm surge
    ! model.
    !---------------------------------------------------
    percentCoriolis = 1.0_SZ
    ! Jie 2014.07
    IF (geof == 1) THEN
      speed = SQRT(((vMax * KT2MS)**2 + vMax * KT2MS * rmx * NM2M * percentCoriolis * corio) *   &
                   (rmx / iDist)**B * EXP(phi * (1.0_SZ - (rmx / iDist)**B)) +                        &
                   (NM2M * iDist * percentCoriolis * corio / 2.0_SZ)**2) -                         &
                  NM2M * iDist * percentCoriolis * corio / 2.0_SZ 
    ELSE 
      speed = SQRT((vMax * KT2MS)**2 * (rmx / iDist)**B * EXP(1.0_SZ - (rmx / iDist)**B) +   &
                   (NM2M * iDist * percentCoriolis * corio / 2.0_SZ)**2) -   &
                  NM2M * iDist * percentCoriolis * corio / 2.0_SZ      
    ENDIF

    ! jgf20111007: Calculate NWS8-like translation speed
    transSpdX = (ABS(speed / (vMax * KT2MS))) * uTrans * KT2MS
    transSpdY = (ABS(speed / (vMax * KT2MS))) * vTrans * KT2MS

    ! now reduce the wind speed to the surface
    speed = speed * windReduction

    u = -speed * COS(DEG2RAD * iAngle)
    v =  speed * SIN(DEG2RAD * iAngle)
    !
    ! Alter wind direction by adding a frictional inflow angle
    CALL Rotate(u, v, FAng(iDist, iRmxTrue), cLat, uf, vf)
    u = uf
    v = vf
    !
    ! jgf20111007: Add in the translation velocity
    u = u + transSpdX
    v = v + transSpdY

    ! convert from 1 minute averaged winds to 10 minute averaged
    ! winds for use in ADCIRC
    u = u * ONE2TEN
    v = v * ONE2TEN

    ! Compute surface pressure from asymmetric hurricane vortex
    IF (geof == 1) THEN
      p = MB2PA * (pc + (pn - pc) * EXP( - phi * (rmx / iDist)**B))
    ELSE
      p = MB2PA * (pc + (pn - pc) * EXP(-(rmx / iDist)**B))
    ENDIF

    ! cut off the vortex field after 401nm
    ! TODO: 401nm should be replaced with something less
    ! arbitrary ... and find a better way to blend this
    !if ( dist > 401.0_SZ ) then
    !u = 0.0_SZ
    !v = 0.0_SZ
    !p = MB2PA * pn
    !endif

  END SUBROUTINE UVPR

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   F A N G
  !----------------------------------------------------------------
  !
  !> Compute a wind angle to parameterize frictional inflow
  !> across isobars.
  !>
  !> On input:
  !>    r           distance from center of storm
  !>    rmx         radius of maximum winds
  !>
  !> On output:
  !>    FAng        frictional inflow angle (degrees)
  !----------------------------------------------------------------
  REAL(SZ) FUNCTION FAng(r, rmx) RESULT(myValOut)

    IMPLICIT NONE

    REAL(SZ), INTENT(IN) :: r
    REAL(SZ), INTENT(IN) :: rmx

    IF ((0.0_SZ <= r) .AND. (r < rmx)) THEN
      myValOut = 10.0_SZ * r / rmx
    ELSE IF ((rmx <= r) .AND. (r < 1.0_SZ * rmx)) THEN
      myValOut = 10.0_SZ + 75.0_SZ * (r / rmx - 1.0_SZ)
    ELSE IF (r >= 1.0_SZ * rmx) THEN
      myValOut = 25.0_SZ
    ELSE
      myValOut = 0.0_SZ
    END IF

  END FUNCTION FAng

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   R O T A T E
  !----------------------------------------------------------------
  !
  !> Rotate a 2D vector (x, y) by an angle.
  !>
  !> On input:
  !>    x           x component of vector
  !>    y           y component of vector
  !>    angle       angle to rotate vector (degrees)
  !>    whichWay    direction of rotation:
  !>                   - = clockwise, + = counter-clockwise
  !>
  !> On output:
  !>    xr          x component of rotated vector
  !>    yr          y component of rotated vector
  !----------------------------------------------------------------
  SUBROUTINE Rotate(x, y, angle, whichWay, xr, yr)

    USE PaHM_Global, ONLY : DEG2RAD

    IMPLICIT NONE

    REAL(SZ), INTENT(IN)  :: x
    REAL(SZ), INTENT(IN)  :: y
    REAL(SZ), INTENT(IN)  :: angle
    REAL(SZ), INTENT(IN)  :: whichWay

    REAL(SZ), INTENT(OUT) :: xr
    REAL(SZ), INTENT(OUT) :: yr

    REAL(SZ)              :: A, cosA, sinA

    A = SIGN(1.0_SZ, whichWay) * DEG2RAD * angle
    cosA = COS(A)
    sinA = SIN(A)

    xr = x * cosA - y * sinA
    yr = x * sinA + y * cosA

  END SUBROUTINE Rotate

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   G E T  L A T E S T  R M A X
  !----------------------------------------------------------------
  REAL(SZ) FUNCTION GetLatestRMax() RESULT(myValOut)

    IMPLICIT NONE

    myValOut = latestRMax

  END FUNCTION GetLatestRMax

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   G E T  L A T E S T  A N G L E
  !----------------------------------------------------------------
  REAL(SZ) FUNCTION GetLatestAngle() RESULT(myValOut)

    IMPLICIT NONE

    myValOut = latestAngle

  END FUNCTION GetLatestAngle

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   V H  W I T H  C O R I  F U L L
  !----------------------------------------------------------------
  !> External function f(x) = 0 for which a root is
  !> sought using Brent's root-finding method.
  !>
  !> On input:
  !>    x       iterative values which converge to root
  !>
  !> On output:
  !>    func    f(x)
  !>
  !> Internal parameters:
  !>    vortex instance variables via accessor functions
  !>
  ! Jie 2013.02 Modified to use the full gradient wind eq.
  !----------------------------------------------------------------
  REAL(SZ) FUNCTION VhWithCoriFull(testRMax) RESULT(myValOut)

    USE PaHM_Global, ONLY : NM2M, KT2MS, MS2KT

    IMPLICIT NONE

    REAL(SZ), INTENT(IN) :: testRMax

    REAL(SZ)             :: thisVR ! the radial wind speed we've been given
    REAL(SZ)             :: vh

    !-------------------------
    ! func(x = rMax) = vh - vr
    !-------------------------
    IF (GetUsequadrantVR() .EQV. .TRUE.) THEN
      thisVR = vrQuadrant(quad)
    ELSE
      thisVR = vr
    END IF

    ! Jie 2013.02
    vh = MS2KT * (SQRT(((vMax * KT2MS)**2 + vMax * KT2MS * testRMax * NM2M * corio) *     &
                       (testRMax / radius(quad))**B *                                        &
                       EXP(phi * (1.0_SZ - (testRMax / radius(quad))**B)) +                   &
                       (NM2M * radius(quad) * corio / 2.0_SZ)**2) -                         &
                  NM2M * radius(quad) * corio / 2.0_SZ)

    myValOut = vh - thisVR

    RETURN

  END FUNCTION VhWithCoriFull

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   V H  W I T H  C O R I
  !----------------------------------------------------------------
  !> External function f(x) = 0 for which a root is
  !> sought using Brent's root-finding method.
  !>
  !> On input:
  !>    x       iterative values which converge to root
  !>
  !> On output:
  !>    func    f(x)
  !>
  !> Internal parameters:
  !>    vortex instance variables via accessor functions
  !----------------------------------------------------------------
  REAL(SZ) FUNCTION VhWithCori(testRMax) RESULT(myValOut)

    USE PaHM_Global, ONLY : NM2M, KT2MS, MS2KT

    IMPLICIT NONE

    REAL(SZ), INTENT(IN) :: testRMax

    REAL(SZ)             :: thisVR ! the radial wind speed we've been given
    REAL(SZ)             :: vh

    !-------------------------
    ! func(x = rMax) = vh - vr
    !-------------------------
    IF (GetUsequadrantVR() .EQV. .TRUE.) THEN
      thisVR = vrQuadrant(quad)
    ELSE
      thisVR = vr
    END IF

    vh = MS2KT * (SQRT((vMax * KT2MS)**2 * (testRMax / radius(quad))**B *                 &
                                              EXP(1.0_SZ - (testRMax / radius(quad))**B) +     &
                       (NM2M * radius(quad) * corio / 2.0_SZ)**2) -                       &
                  NM2M * radius(quad) * corio / 2.0_SZ)

    myValOut = vh - thisVR

    RETURN

  END FUNCTION VhWithCori

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   V H  N O  C O R I
  !----------------------------------------------------------------
  REAL(SZ) FUNCTION VhNoCori(testRMax) RESULT(myValOut)

    USE PaHM_Global, ONLY : KT2MS, MS2KT

    IMPLICIT NONE

    REAL(SZ), INTENT(IN) :: testRMax

    REAL(SZ) :: thisVR ! the radial wind speed we've been given

    IF (GetUsequadrantVR() .EQV. .TRUE.) THEN
      thisVR = vrQuadrant(quad)
    ELSE
      thisVR = vr
    END IF

    myValOut = ABS(MS2KT * SQRT((vMax * KT2MS)**2 * (testRMax / radius(quad))**B *     &
                                EXP(1 - (testRMax / radius(quad))**B))) - thisVR

    RETURN

  END FUNCTION VhNoCori

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   F I N D  R O O T
  !----------------------------------------------------------------
  !> Use brute-force marching to find a root the interval [x1,x2].
  !>
  !> On input:
  !>    func        function f(x)=0 for which root is sough
  !>    x1          left  side of interval
  !>    x2          right side of interval
  !>    dx          x increment for march
  !>
  !> On output:
  !>    a           left  side of interval that brackets the root
  !>    b           right side of interval that brackets the root
  !>    FindRoot    root returned
  !----------------------------------------------------------------
  REAL(SZ) FUNCTION FindRoot(func, x1, x2, dx, a, b) RESULT(myRoot)
!PV Need to check for the x2 variable is not used anywhere next
    IMPLICIT NONE

    REAL(SZ), EXTERNAL    :: func
    REAL(SZ), INTENT(IN)  :: x1, x2          ! Search interval [x1,x2]
    REAL(SZ), INTENT(IN)  :: dx              ! Marching increment
    REAL(SZ), INTENT(OUT) :: a, b            ! x values that bracket root

    INTEGER , PARAMETER   :: ITERMAX = 400   ! Max # of iterations
    INTEGER               :: iter            ! iteration counter
    REAL(SZ)              :: fa, fb          ! function values f(x)

    ! Initialize left side of interval
    a  = x1
    fa = func(a)

    ! March along interval until root is found
    ! or solution diverges.
    myRoot = a
    DO iter = 1, ITERMAX
      b  = x1 + iter * dx
      fb = func(b)

      ! Check progress
      IF ((fa * fb < 0.0_SZ) .OR. (ABS(fb) > ABS(fa))) THEN
        ! Assign root
        IF (ABS(fb) > ABS(fa)) THEN
          myRoot = a
        ELSE
          myRoot = b
        END IF

        EXIT
      END IF

      ! Move right search interval values to left side
      ! for next iteration.
      a  = b
      fa = fb
    END DO

    IF (iter >= ITERMAX) THEN
      PRINT *, "FUNCTION FindRoot: exceeded max # of iterations"
      myRoot = -99999.0
    END IF

    RETURN

  END FUNCTION FindRoot

!================================================================================

END MODULE PaHM_Vortex
