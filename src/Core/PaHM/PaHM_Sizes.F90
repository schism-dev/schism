!----------------------------------------------------------------
!               M O D U L E   S I Z E S
!----------------------------------------------------------------
!> @file sizes.F90
!>
!> @brief
!>    Contains the definitions of various number types and utilities used in PaHM.
!>
!> @details
!>   
!>
!> @author Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>
!----------------------------------------------------------------

MODULE PaHM_Sizes

  IMPLICIT NONE

  !-----------------------------------------------------------------------
  ! I N T E R F A C E S
  !-----------------------------------------------------------------------
  INTERFACE CompareReals
    MODULE PROCEDURE CompareSingleReals
    MODULE PROCEDURE CompareDoubleReals
  END INTERFACE CompareReals

  INTERFACE FixNearWholeReal
    MODULE PROCEDURE FixNearWholeSingleReal
    MODULE PROCEDURE FixNearWholeDoubleReal
  END INTERFACE FixNearWholeReal
  !-----------------------------------------------------------------------

  ! SP = single precision, HP = high (double) precision
  INTEGER, PARAMETER :: SP = SELECTED_REAL_KIND(6,   37)   !  6  digits, range \([10^{-37}  , 10^{+37}   - 1]\), 32 bits
  INTEGER, PARAMETER :: HP = SELECTED_REAL_KIND(15, 307)   ! 15  digits, range \([10^{-307} , 10^{+307}  - 1]\), 64 bits

  ! Precision of integers:
  INTEGER, PARAMETER :: INT16 = SELECTED_INT_KIND(38)      ! Range \([-2^{127},+2^{127} - 1]\), 39 digits plus sign; 128 bits
  INTEGER, PARAMETER ::  INT8 = SELECTED_INT_KIND(18)      ! Range \([-2^{63},+2^{63} - 1]\),   19 digits plus sign;  64 bits
  INTEGER, PARAMETER ::  INT4 = SELECTED_INT_KIND( 9)      ! Range \([-2^{31},+2^{31} - 1]\),   10 digits plus sign;  32 bits
  INTEGER, PARAMETER ::  INT2 = SELECTED_INT_KIND( 4)      ! Range \([-2^{15},+2^{15} - 1]\),    5 digits plus sign;  16 bits
  INTEGER, PARAMETER ::  INT1 = SELECTED_INT_KIND( 2)      ! Range \([-2^{7} ,+2^{7}  - 1]\),    3 digits plus sign;   8 bits
  INTEGER, PARAMETER :: LONG  = INT8
  INTEGER, PARAMETER :: LLONG = INT16

  INTEGER,PARAMETER :: WP = HP   ! default real kind (for csv_module)
  INTEGER,PARAMETER :: IP = INT8 ! default integer kind (for csv_module)

  ! By default we perform all calculations in double precision
  ! SET NUMBER OF BYTES "SZ" IN REAL(SZ) DECLARATIONS
  ! SET "NBYTE" FOR PROCESSING INPUT DATA RECORD LENGTH
#ifdef REAL4
  INTEGER, PARAMETER :: SZ = SP
  INTEGER, PARAMETER :: NBYTE = 4
#else
  INTEGER, PARAMETER :: SZ = HP
  INTEGER, PARAMETER :: NBYTE = 8
#endif

  ! Used to initialize the mesh arrays and in NetCDF output files for missing values.
  ! Also used to initialize some input variables to check if these variables
  ! were supplied user defined values.
  REAL(SZ), PARAMETER :: RMISSV = -999999.0_SZ
  INTEGER, PARAMETER  :: IMISSV = -999999

  CHARACTER(LEN=1), PARAMETER :: BLANK = ' '
  
  ! Filename length (considers the presence of the full path in the filename)
  INTEGER, PARAMETER :: FNAMELEN = 1024


  CONTAINS


  !----------------------------------------------------------------
  ! F U N C T I O N   CompareDoubleReals
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Compares two double precision numbers.
  !>
  !> @details
  !>   Allow users to define the value of eps. If not, eps equals to the default machine eps.
  !>
  !> @param[in]
  !>   rVal1   The first value (double precision number) in the comparison
  !> @param[in]
  !>   rVal2   The second value (double precision number) in the comparison
  !> @param[in]
  !>   eps     The tolerance (optional) for the comparison
  !>
  !> @return myValOut
  !> @verbatim
  !>   -1 (if rVal1 < rVal2)
  !>    0 (if rVal1 = rVal2)
  !>   +1 (if rVal1 > rVal2)
  !> @endverbatim
  !>
  !> @note The code was adopted from the D-Flow FM source (...src/precision_basics.F90)
  !----------------------------------------------------------------
  INTEGER FUNCTION CompareDoubleReals(rVal1, rVal2, eps) RESULT(myValOut)

    IMPLICIT NONE

    ! Global variables
    REAL(HP), INTENT(IN)           :: rVal1, rVal2
    REAL(HP), OPTIONAL, INTENT(IN) :: eps

    ! Local variables
    REAL(HP)                       :: epsSys, epsUsr, value


    epsSys = 2.0_HP * EPSILON(rVal1)

    IF (PRESENT(eps)) THEN
      epsUsr = ABS(eps)
    ELSE
      epsUsr = epsSys
    ENDIF

    IF ((ABS(rVal1) < 1.0_HP) .OR. (ABS(rVal2) < 1.0_HP)) THEN
      value = rVal1 - rVal2
    ELSE
      value = (rVal1 - rVal2) / MAX(rVal1, rVal2)
      IF (ABS(value) < 1.0_HP) value = rVal1 - rVal2
    END IF

    IF (ABS(value) < epsUsr) THEN
      myValOut = 0
    ELSE IF (rVal1 < rVal2) THEN
      myValOut = -1
    ELSE
      myValOut = 1
    END IF

    RETURN

  END FUNCTION CompareDoubleReals

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   C O M P A R E  S I N G L E  R E A L S
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Compares two single precision numbers.
  !>
  !> @details
  !>   Allow users to define the value of eps. If not, eps equals to the default machine eps.
  !>
  !> @param[in]
  !>   rVal1   The first value (single precision number) in the comparison
  !> @param[in]
  !>   rVal2   The second value (single precision number) in the comparison
  !> @param[in]
  !>   eps     The tolerance (optional) for the comparison
  !>
  !> @return myValOut
  !> @verbatim
  !>   -1 (if rVal1 < rVal2)
  !>    0 (if rVal1 = rVal2)
  !>   +1 (if rVal1 > rVal2)
  !> @endverbatim
  !>
  !> @note The code was adopted from the D-Flow FM source (...src/precision_basics.F90)
  !----------------------------------------------------------------
  INTEGER FUNCTION CompareSingleReals(rVal1, rVal2, eps) RESULT(myValOut)

    IMPLICIT NONE

    ! Global variables
    REAL(SP), INTENT(IN)           :: rVal1, rVal2
    REAL(SP), OPTIONAL, INTENT(IN) :: eps

    ! Local variables
    REAL(SP)                       :: epsSys, epsUsr, value


    epsSys = 2.0_SP * EPSILON(rVal1)

    IF (PRESENT(eps)) THEN
      epsUsr = ABS(eps)
    ELSE
      epsUsr = epsSys
    ENDIF

    IF ((ABS(rVal1) < 1.0_SP) .OR. (ABS(rVal2) < 1.0_SP)) THEN
      value = rVal1 - rVal2
    ELSE
      value = (rVal1 - rVal2) / MAX(rVal1, rVal2)
      IF (ABS(value) < 1.0_SP) value = rVal1 - rVal2
    END IF

    IF (ABS(value) < epsUsr) THEN
      myValOut = 0
    ELSE IF (rVal1 < rVal2) THEN
      myValOut = -1
    ELSE
      myValOut = 1
    END IF

    RETURN

  END FUNCTION CompareSingleReals

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   F I X  N E A R  W H O L E  D O U B L E  R E A L
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Rounds a double precision real number to its nearest whole number.
  !>
  !> @details
  !>   Rounds a double precision real number to its nearest whole number.
  !>   If the real number is very close (within a tolerance) to its nearest whole number
  !>   then it is set equal to its nearest whole number.
  !>   Allow users to define the value of the tolerance "eps". If not, then eps equals
  !>   to the default machine eps.
  !>
  !> @param[in]
  !>   rVal   The real number value (double precision) in the comparison
  !> @param[in]
  !>   eps    The tolerance (optional) for the comparison
  !>
  !> @return myValOut : Either **rVal** or its nearest integer **iVar** converted to double
  !> @verbatim
  !>   rVal (if abs(rVal - iVal) >  eps
  !>   iVal (if abs(rVal - iVal) <= eps
  !> @endverbatim
  !>
  !----------------------------------------------------------------
  REAL(HP) FUNCTION FixNearWholeDoubleReal(rVal, eps) RESULT(myValOut)

    IMPLICIT NONE

    ! Global Variables
    REAL(HP), INTENT(IN)           :: rVal
    REAL(HP), OPTIONAL, INTENT(IN) :: eps

    ! Local Variables
    REAL(HP)                       :: epsSys, epsUsr, value


    epsSys = 2.0_HP * EPSILON(rVal)

    IF (PRESENT(eps)) THEN
      epsUsr = ABS(eps)
    ELSE
      epsUsr = epsSys
    ENDIF

    myValOut = rVal
    value    = ANINT(myValOut)
    IF (CompareReals(myValOut, value, epsUsr) == 0) myValOut = value

    RETURN

  END FUNCTION FixNearWholeDoubleReal

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   F I X  N E A R  W H O L E  S I N G L E  R E A L
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Rounds a single precision real number to its nearest whole number.
  !>
  !> @details
  !>   Rounds a single precision real number to its nearest whole number.
  !>   If the real number is very close (within a tolerance) to its nearest whole number
  !>   then it is set equal to its nearest whole number.
  !>   Allow users to define the value of the tolerance "eps". If not, then eps equals
  !>   to the default machine eps.
  !>
  !> @param[in]
  !>   rVal   The real number value (single precision) in the comparison
  !> @param[in]
  !>   eps    The tolerance (optional) for the comparison
  !>
  !> @return myValOut : Either **rVal** or its nearest integer **iVar** converted to real
  !> @verbatim
  !>   rVal (if abs(rVal - iVal) >  eps
  !>   iVal (if abs(rVal - iVal) <= eps
  !> @endverbatim
  !>
  !----------------------------------------------------------------
  REAL(SP) FUNCTION FixNearWholeSingleReal(rVal, eps) RESULT(myValOut)

    IMPLICIT NONE

    ! Global Variables
    REAL(SP), INTENT(IN)           :: rVal
    REAL(SP), OPTIONAL, INTENT(IN) :: eps

    ! Local Variables
    REAL(SP)                       :: epsSys, epsUsr, value


    epsSys = 2.0_SP * EPSILON(rVal)

    IF (PRESENT(eps)) THEN
      epsUsr = ABS(eps)
    ELSE
      epsUsr = epsSys
    ENDIF

    myValOut = rVal
    value    = ANINT(myValOut)
    IF (CompareReals(myValOut, value, epsUsr) == 0) myValOut = value

    RETURN

  END FUNCTION FixNearWholeSingleReal

!================================================================================

END MODULE PaHM_Sizes
