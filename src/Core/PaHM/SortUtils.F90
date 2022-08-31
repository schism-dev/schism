!----------------------------------------------------------------
!               M O D U L E   U T I L I T I E S
!----------------------------------------------------------------
!> @file sortutils.F90
!>
!>
!> @brief
!>   
!>
!> @details
!>   
!>
!> @author Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>
!----------------------------------------------------------------

MODULE SortUtils

  USE PaHM_Sizes
  USE PaHM_Messages

  !-----------------------------------------------------------------------
  ! I N T E R F A C E S
  !-----------------------------------------------------------------------
  INTERFACE Indexx
    MODULE PROCEDURE IndexxInt
    MODULE PROCEDURE IndexxInt8
    MODULE PROCEDURE IndexxString
    MODULE PROCEDURE IndexxSingle
    MODULE PROCEDURE IndexxDouble
  END INTERFACE Indexx

  INTERFACE Arth
    MODULE PROCEDURE ArthInt
    MODULE PROCEDURE ArthSingle
    MODULE PROCEDURE ArthDouble
  END INTERFACE Arth

  INTERFACE ArrayCopy
    MODULE PROCEDURE ArrayCopyInt
    MODULE PROCEDURE ArrayCopySingle
    MODULE PROCEDURE ArrayCopyDouble
  END INTERFACE ArrayCopy

  INTERFACE ArrayEqual
    MODULE PROCEDURE ArrayEqualInt
    MODULE PROCEDURE ArrayEqualSingle
    MODULE PROCEDURE ArrayEqualDouble
  END INTERFACE ArrayEqual

  INTERFACE Swap
    MODULE PROCEDURE SwapInt
    MODULE PROCEDURE SwapSingle
    MODULE PROCEDURE SwapDouble
    MODULE PROCEDURE SwapIntVec
    MODULE PROCEDURE SwapSingleVec
    MODULE PROCEDURE SwapDoubleVec
  END INTERFACE Swap
  !-----------------------------------------------------------------------


  CONTAINS


  !----------------------------------------------------------------
  ! S U B R O U T I N E   I N D E X X  I N T
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Indexes a 1D integer array in ascending order.
  !>
  !> @details
  !>   Indexes the 1D array arr1D, i.e., outputs the array index of length N such that arr1D(idx1D(j ))
  !>   is in ascending order for j = 1, 2, . . . , N. The input quantity arr1D is not changed.
  !>
  !> @param[in]
  !>   arr1D    The array to be indexed (integer)
  !> @param[out]
  !>   idx1D    The array of "indexed" indexes of arr1D (output)
  !> @param[out]
  !>   status   The error status, no error: status = 0 (output)
  !>
  !> @note Adopted from Numerical Recipes for Fortran 90
  !>
  !----------------------------------------------------------------
  SUBROUTINE IndexxInt(arr1D, idx1D, status)

    IMPLICIT NONE

    ! Global variables
    INTEGER, DIMENSION(:), INTENT(IN)  :: arr1D
    INTEGER, DIMENSION(:), INTENT(OUT) :: idx1D
    INTEGER, OPTIONAL, INTENT(OUT)     :: status

    ! Local variables
    INTEGER, PARAMETER                 :: NN = 15, NSTACK = 50
    INTEGER                            :: a
    INTEGER                            :: nARR, nIDX, tmpIDX
    INTEGER                            :: k, i, j, l, r
    INTEGER                            :: ist, stack(NSTACK)
    CHARACTER(LEN=64)                  :: tmpStr1, tmpStr2


    CALL SetMessageSource("IndexxInt")

    IF (PRESENT(status)) status = 0

    nARR = SIZE(arr1D, 1)
    nIDX = SIZE(idx1D, 1)

    IF (nARR /= nIDX) THEN
      WRITE(tmpStr1, '(a, i0)') 'nARR = ', nARR
      WRITE(tmpStr2, '(a, i0)') 'nIDX = ', nIDX
      WRITE(scratchMessage, '(a)') 'The size of the 1D arrays arr1D and idx1D is not the same: ' // &
                                   TRIM(ADJUSTL(tmpStr1)) // ', ' // TRIM(ADJUSTL(tmpStr2))
      
      CALL AllMessage(ERROR, scratchMessage)
      CALL UnsetMessageSource()

      IF (PRESENT(status)) status = 1

      RETURN
    END IF

    idx1D = Arth(1, 1, nARR)

    ist = 0
    l   = 1
    r   = nARR

    DO
      IF (r - l < NN) THEN
        DO j = l + 1, r
          tmpIDX = idx1D(j)
          a = arr1D(tmpIDX)
          DO i = j - 1, l, -1
            IF (arr1D(idx1D(i)) <= a) EXIT
            idx1D(i + 1) = idx1D(i)
          END DO
          idx1D(i + 1) = tmpIDX
        END DO

        IF (ist == 0) THEN
          CALL UnsetMessageSource()

          RETURN
        END IF

        r   = stack(ist)
        l   = stack(ist - 1)
        ist = ist - 2
      ELSE
        k = (l + r) / 2

        CALL Swap(idx1D(k), idx1D(l + 1))
        CALL IcompXchg(idx1D(l), idx1D(r))
        CALL IcompXchg(idx1D(l + 1), idx1D(r))
        CALL IcompXchg(idx1D(l), idx1D(l + 1))

        i = l + 1
        j = r
        tmpIDX = idx1D(l + 1)
        a = arr1D(tmpIDX)

        DO
          DO
            i = i + 1
            IF (arr1D(idx1D(i)) > a) EXIT
          END DO

          DO
            j = j - 1
            IF (arr1D(idx1D(j)) < a) EXIT
          END DO

          IF (j < i) EXIT
          CALL Swap(idx1D(i), idx1D(j))
        END DO

        idx1D(l + 1) = idx1D(j)
        idx1D(j) = tmpIDX
        ist = ist + 2

        IF (ist > NSTACK) THEN
          WRITE(tmpStr1, '(a, i0)') 'NSTACK = ', NSTACK
          WRITE(scratchMessage, '(a)') 'The value of the NSTACK parameter is too small: ' // &
                                       TRIM(ADJUSTL(tmpStr1))

          CALL LogMessage(ERROR, scratchMessage)
          CALL UnsetMessageSource()

          IF (PRESENT(status)) status = 2

          RETURN

        END IF

        IF (r - i + 1 >= j - l) THEN
          stack(ist) = r
          stack(ist - 1) = i
          r = j - 1
        ELSE
          stack(ist) = j - 1
          stack(ist - 1) = l
          l = i
        END IF
      END IF
    END DO

    CALL UnSetMessageSource()


    CONTAINS

    SUBROUTINE IcompXchg(i, j)

      IMPLICIT NONE

      ! Global variables
      INTEGER, INTENT(INOUT) :: i, j

      ! Local variables
      INTEGER :: swp

      IF (arr1D(j) < arr1D(i)) THEN
        swp = i
        i   = j
        j   = swp
      END IF

    END SUBROUTINE IcompXchg

  END SUBROUTINE IndexxInt

  !================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   I N D E X X  I N T  8
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Indexes a 1D 32-bit integer array in ascending order.
  !>
  !> @details
  !>   Indexes the 1D array arr1D, i.e., outputs the array index of length N such that arr1D(idx1D(j ))
  !>   is in ascending order for j = 1, 2, . . . , N. The input quantity arr1D is not changed.
  !>
  !> @param[in]
  !>   arr1D    The array to be indexed (integer)
  !> @param[out]
  !>   idx1D    The array of "indexed" indexes of arr1D (output)
  !> @param[out]
  !>   status   The error status, no error: status = 0 (output)
  !>
  !> @note Adopted from Numerical Recipes for Fortran 90
  !>
  !----------------------------------------------------------------
  SUBROUTINE IndexxInt8(arr1D, idx1D, status)

    IMPLICIT NONE

    ! Global variables
    INTEGER(INT8), DIMENSION(:), INTENT(IN)  :: arr1D
    INTEGER, DIMENSION(:), INTENT(OUT)       :: idx1D
    INTEGER, OPTIONAL, INTENT(OUT)           :: status

    ! Local variables
    INTEGER, PARAMETER                       :: NN = 15, NSTACK = 50
    INTEGER(INT8)                            :: a
    INTEGER                                  :: nARR, nIDX, tmpIDX
    INTEGER                                  :: k, i, j, l, r
    INTEGER                                  :: ist, stack(NSTACK)
    CHARACTER(LEN=64)                        :: tmpStr1, tmpStr2


    CALL SetMessageSource("IndexxInt8")

    IF (PRESENT(status)) status = 0

    nARR = SIZE(arr1D, 1)
    nIDX = SIZE(idx1D, 1)

    IF (nARR /= nIDX) THEN
      WRITE(tmpStr1, '(a, i0)') 'nARR = ', nARR
      WRITE(tmpStr2, '(a, i0)') 'nIDX = ', nIDX
      WRITE(scratchMessage, '(a)') 'The size of the 1D arrays arr1D and idx1D is not the same: ' // &
                                   TRIM(ADJUSTL(tmpStr1)) // ', ' // TRIM(ADJUSTL(tmpStr2))
      
      CALL AllMessage(ERROR, scratchMessage)
      CALL UnsetMessageSource()

      IF (PRESENT(status)) status = 1

      RETURN
    END IF

    idx1D = Arth(1, 1, nARR)

    ist = 0
    l   = 1
    r   = nARR

    DO
      IF (r - l < NN) THEN
        DO j = l + 1, r
          tmpIDX = idx1D(j)
          a = arr1D(tmpIDX)
          DO i = j - 1, l, -1
            IF (arr1D(idx1D(i)) <= a) EXIT
            idx1D(i + 1) = idx1D(i)
          END DO
          idx1D(i + 1) = tmpIDX
        END DO

        IF (ist == 0) THEN
          CALL UnsetMessageSource()

          RETURN
        END IF

        r   = stack(ist)
        l   = stack(ist - 1)
        ist = ist - 2
      ELSE
        k = (l + r) / 2

        CALL Swap(idx1D(k), idx1D(l + 1))
        CALL IcompXchg(idx1D(l), idx1D(r))
        CALL IcompXchg(idx1D(l + 1), idx1D(r))
        CALL IcompXchg(idx1D(l), idx1D(l + 1))

        i = l + 1
        j = r
        tmpIDX = idx1D(l + 1)
        a = arr1D(tmpIDX)

        DO
          DO
            i = i + 1
            IF (arr1D(idx1D(i)) > a) EXIT
          END DO

          DO
            j = j - 1
            IF (arr1D(idx1D(j)) < a) EXIT
          END DO

          IF (j < i) EXIT
          CALL Swap(idx1D(i), idx1D(j))
        END DO

        idx1D(l + 1) = idx1D(j)
        idx1D(j) = tmpIDX
        ist = ist + 2

        IF (ist > NSTACK) THEN
          WRITE(tmpStr1, '(a, i0)') 'NSTACK = ', NSTACK
          WRITE(scratchMessage, '(a)') 'The value of the NSTACK parameter is too small: ' // &
                                       TRIM(ADJUSTL(tmpStr1))

          CALL LogMessage(ERROR, scratchMessage)
          CALL UnsetMessageSource()

          IF (PRESENT(status)) status = 2

          RETURN

        END IF

        IF (r - i + 1 >= j - l) THEN
          stack(ist) = r
          stack(ist - 1) = i
          r = j - 1
        ELSE
          stack(ist) = j - 1
          stack(ist - 1) = l
          l = i
        END IF
      END IF
    END DO

    CALL UnSetMessageSource()


    CONTAINS

    SUBROUTINE IcompXchg(i, j)

      IMPLICIT NONE

      ! Global variables
      INTEGER, INTENT(INOUT) :: i, j

      ! Local variables
      INTEGER :: swp

      IF (arr1D(j) < arr1D(i)) THEN
        swp = i
        i   = j
        j   = swp
      END IF

    END SUBROUTINE IcompXchg

  END SUBROUTINE IndexxInt8

  !================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   I N D E X X  S T R I N G
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Indexes a 1D string array in ascending order.
  !>
  !> @details
  !>   Indexes the 1D array arr1D, i.e., outputs the array index of length N such that arr1D(idx1D(j ))
  !>   is in ascending order for j = 1, 2, . . . , N. The input quantity arr1D is not changed.
  !>   Modified version of IndexxInt to account for string comparisons
  !>
  !> @param[in]
  !>   arr1D      The array to be indexed (string)
  !> @param[out]
  !>   idx1D      The array of "indexed" indexes of arr1D (output)
  !> @param[out]
  !>   status     The error status, no error: status = 0 (output)
  !> @param[in]
  !>   caseSens   Logical flag to request case sensitive sort
  !>
  !----------------------------------------------------------------
  SUBROUTINE IndexxString(arr1D, idx1D, status, caseSens)

    IMPLICIT NONE

    ! Global variables
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: arr1D
    LOGICAL, OPTIONAL, INTENT(IN)              :: caseSens
    INTEGER, DIMENSION(:), INTENT(OUT)         :: idx1D
    INTEGER, OPTIONAL, INTENT(OUT)             :: status

    ! Local variables
    INTEGER, PARAMETER                         :: NN = 15, NSTACK = 50
    CHARACTER(LEN=LEN(arr1D(1)))               :: a
    INTEGER                                    :: nARR, nIDX, tmpIDX
    INTEGER                                    :: k, i, j, l, r
    INTEGER                                    :: ist, stack(NSTACK)
    CHARACTER(LEN=64)                          :: tmpStr1, tmpStr2
    LOGICAL                                    :: sFlag


    CALL SetMessageSource("IndexxString")

    sFlag = .TRUE.
    IF (PRESENT(caseSens)) sFlag = caseSens

    IF (PRESENT(status)) status   = 0

    nARR = SIZE(arr1D, 1)
    nIDX = SIZE(idx1D, 1)

    IF (nARR /= nIDX) THEN
      WRITE(tmpStr1, '(a, i0)') 'nARR = ', nARR
      WRITE(tmpStr2, '(a, i0)') 'nIDX = ', nIDX
      WRITE(scratchMessage, '(a)') 'The size of the 1D arrays arr1D and idx1D is not the same: ' // &
                                   TRIM(ADJUSTL(tmpStr1)) // ', ' // TRIM(ADJUSTL(tmpStr2))
      
      CALL AllMessage(ERROR, scratchMessage)
      CALL UnsetMessageSource()

      IF (PRESENT(status)) status = 1

      RETURN
    END IF

    idx1D = Arth(1, 1, nARR)

    ist = 0
    l   = 1
    r   = nARR

    DO
      IF (r - l < NN) THEN
        DO j = l + 1, r
          tmpIDX = idx1D(j)
          a = arr1D(tmpIDX)
          DO i = j - 1, l, -1
            IF (StringLexComp(arr1D(idx1D(i)), a, sFlag) <= 0) EXIT
            idx1D(i + 1) = idx1D(i)
          END DO
          idx1D(i + 1) = tmpIDX
        END DO

        IF (ist == 0) THEN
          CALL UnsetMessageSource()

          RETURN
        END IF

        r   = stack(ist)
        l   = stack(ist - 1)
        ist = ist - 2
      ELSE
        k = (l + r) / 2

        CALL Swap(idx1D(k), idx1D(l + 1))
        CALL IcompXchg(idx1D(l), idx1D(r))
        CALL IcompXchg(idx1D(l + 1), idx1D(r))
        CALL IcompXchg(idx1D(l), idx1D(l + 1))
          
        i = l + 1
        j = r
        tmpIDX = idx1D(l + 1)
        a = arr1D(tmpIDX)

        DO
          DO
            i = i + 1
            IF (StringLexComp(arr1D(idx1D(i)), a, sFlag) > 0) EXIT
          END DO

          DO
            j = j - 1
            IF (StringLexComp(arr1D(idx1D(j)), a, sFlag) < 0) EXIT
          END DO

          IF (j < i) EXIT
          CALL Swap(idx1D(i), idx1D(j))
        END DO

        idx1D(l + 1) = idx1D(j)
        idx1D(j) = tmpIDX
        ist = ist + 2

        IF (ist > NSTACK) THEN
          WRITE(tmpStr1, '(a, i0)') 'NSTACK = ', NSTACK
          WRITE(scratchMessage, '(a)') 'The value of the NSTACK parameter is too small: ' // &
                                       TRIM(ADJUSTL(tmpStr1))

          CALL LogMessage(ERROR, scratchMessage)
          CALL UnsetMessageSource()

          IF (PRESENT(status)) status = 2

          RETURN

        END IF

        IF (r - i + 1 >= j - l) THEN
          stack(ist) = r
          stack(ist - 1) = i
          r = j - 1
        ELSE
          stack(ist) = j - 1
          stack(ist - 1) = l
          l = i
        END IF
      END IF
    END DO

    CALL UnsetMessageSource()


    CONTAINS

    SUBROUTINE IcompXchg(i, j)

      IMPLICIT NONE

      ! Global variables
      INTEGER, INTENT(INOUT) :: i, j

      ! Local variables
      INTEGER :: swp

      IF (StringLexComp(arr1D(j), arr1D(i), sFlag) < 0) THEN
        swp = i
        i   = j
        j   = swp
      END IF

    END SUBROUTINE IcompXchg

  END SUBROUTINE IndexxString

  !================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   I N D E X X  S I N G L E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Indexes a 1D single precision array in ascending order.
  !>
  !> @details
  !>   Indexes the 1D array arr1D, i.e., outputs the array index of length N such that arr1D(idx1D(j ))
  !>   is in ascending order for j = 1, 2, . . . , N. The input quantity arr1D is not changed.
  !>
  !> @param[in]
  !>   arr1D      The array to be indexed (single precision)
  !> @param[out]
  !>   idx1D      The array of "indexed" indexes of arr1D (output)
  !> @param[out]
  !>   status     The error status, no error: status = 0 (output)
  !>
  !> @note Adopted from Numerical Recipes for Fortran 90
  !>
  !----------------------------------------------------------------
  SUBROUTINE IndexxSingle(arr1D, idx1D, status)

    IMPLICIT NONE

    ! Global variables
    REAL(SP), DIMENSION(:), INTENT(IN) :: arr1D
    INTEGER, DIMENSION(:), INTENT(OUT) :: idx1D
    INTEGER, OPTIONAL, INTENT(OUT)     :: status

    ! Local variables
    INTEGER, PARAMETER                 :: NN = 15, NSTACK = 50
    REAL(SP)                           :: a
    INTEGER                            :: nARR, nIDX, tmpIDX
    INTEGER                            :: k, i, j, l, r
    INTEGER                            :: ist, stack(NSTACK)
    CHARACTER(LEN=64)                  :: tmpStr1, tmpStr2


    CALL SetMessageSource("IndexxSingle")

    IF (PRESENT(status)) status = 0

    nARR = SIZE(arr1D, 1)
    nIDX = SIZE(idx1D, 1)

    IF (nARR /= nIDX) THEN
      WRITE(tmpStr1, '(a, i0)') 'nARR = ', nARR
      WRITE(tmpStr2, '(a, i0)') 'nIDX = ', nIDX
      WRITE(scratchMessage, '(a)') 'The size of the 1D arrays arr1D and idx1D is not the same: ' // &
                                   TRIM(ADJUSTL(tmpStr1)) // ', ' // TRIM(ADJUSTL(tmpStr2))
      
      CALL LogMessage(ERROR, scratchMessage)
      CALL UnsetMessageSource()

      IF (PRESENT(status)) status = 1

      RETURN
    END IF

    idx1D = Arth(1, 1, nARR)

    ist = 0
    l   = 1
    r   = nARR

    DO
      IF (r - l < NN) THEN
        DO j = l + 1, r
          tmpIDX = idx1D(j)
          a = arr1D(tmpIDX)
          DO i = j - 1, l, -1
            IF (arr1D(idx1D(i)) <= a) EXIT
            idx1D(i + 1) = idx1D(i)
          END DO
          idx1D(i + 1) = tmpIDX
        END DO

        IF (ist == 0) THEN
          CALL UnsetMessageSource()

          RETURN
        END IF

        r   = stack(ist)
        l   = stack(ist - 1)
        ist = ist - 2
      ELSE
        k = (l + r) / 2

        CALL Swap(idx1D(k), idx1D(l + 1))
        CALL IcompXchg(idx1D(l), idx1D(r))
        CALL IcompXchg(idx1D(l + 1), idx1D(r))
        CALL IcompXchg(idx1D(l), idx1D(l + 1))

        i = l + 1
        j = r
        tmpIDX = idx1D(l + 1)
        a = arr1D(tmpIDX)

        DO
          DO
            i = i + 1
            IF (arr1D(idx1D(i)) > a) EXIT
          END DO

          DO
            j = j - 1
            IF (arr1D(idx1D(j)) < a) EXIT
          END DO

          IF (j < i) EXIT
          CALL Swap(idx1D(i), idx1D(j))
        END DO

        idx1D(l + 1) = idx1D(j)
        idx1D(j) = tmpIDX
        ist = ist + 2

        IF (ist > NSTACK) THEN
          WRITE(tmpStr1, '(a, i0)') 'NSTACK = ', NSTACK
          WRITE(scratchMessage, '(a)') 'The value of the NSTACK parameter is too small: ' // &
                                       TRIM(ADJUSTL(tmpStr1))

          CALL LogMessage(ERROR, scratchMessage)
          CALL UnsetMessageSource()

          IF (PRESENT(status)) status = 2

          RETURN

        END IF

        IF (r - i + 1 >= j - l) THEN
          stack(ist) = r
          stack(ist - 1) = i
          r = j - 1
        ELSE
          stack(ist) = j - 1
          stack(ist - 1) = l
          l = i
        END IF
      END IF
    END DO

    CALL UnsetMessageSource()


    CONTAINS

    SUBROUTINE IcompXchg(i, j)

      IMPLICIT NONE

      ! Global variables
      INTEGER, INTENT(INOUT) :: i, j

      ! Local variables
      INTEGER :: swp

      IF (arr1D(j) < arr1D(i)) THEN
        swp = i
        i   = j
        j   = swp
      END IF

    END SUBROUTINE IcompXchg

  END SUBROUTINE IndexxSingle

  !================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   I N D E X X  D O U B L E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Indexes a 1D double precision array in ascending order.
  !>
  !> @details
  !>   Indexes the 1D array arr1D, i.e., outputs the array index of length N such that arr1D(idx1D(j ))
  !>   is in ascending order for j = 1, 2, . . . , N. The input quantity arr1D is not changed.
  !>
  !> @param[in]
  !>   arr1D      The array to be indexed (double precision)
  !> @param[out]
  !>   idx1D      The array of "indexed" indexes of arr1D (output)
  !> @param[out]
  !>   status     The error status, no error: status = 0 (output)
  !>
  !> @note Adopted from Numerical Recipes for Fortran 90
  !>
  !----------------------------------------------------------------
  SUBROUTINE IndexxDouble(arr1D, idx1D, status)

    IMPLICIT NONE

    ! Global variables
    REAL(HP), DIMENSION(:), INTENT(IN) :: arr1D
    INTEGER, DIMENSION(:), INTENT(OUT) :: idx1D
    INTEGER, OPTIONAL, INTENT(OUT)     :: status

    ! Local variables
    INTEGER, PARAMETER                 :: NN = 15, NSTACK = 50
    REAL(HP)                           :: a
    INTEGER                            :: nARR, nIDX, tmpIDX
    INTEGER                            :: k, i, j, l, r
    INTEGER                            :: ist, stack(NSTACK)
    CHARACTER(LEN=64)                  :: tmpStr1, tmpStr2


    CALL SetMessageSource("IndexxDouble")

    IF (PRESENT(status)) status = 0

    nARR = SIZE(arr1D, 1)
    nIDX = SIZE(idx1D, 1)

    IF (nARR /= nIDX) THEN
      WRITE(tmpStr1, '(a, i0)') 'nARR = ', nARR
      WRITE(tmpStr2, '(a, i0)') 'nIDX = ', nIDX
      WRITE(scratchMessage, '(a)') 'The size of the 1D arrays arr1D and idx1D is not the same: ' // &
                                   TRIM(ADJUSTL(tmpStr1)) // ', ' // TRIM(ADJUSTL(tmpStr2))
      
      CALL LogMessage(ERROR, scratchMessage)
      CALL UnsetMessageSource()

      IF (PRESENT(status)) status = 1

      RETURN
    END IF

    idx1D = Arth(1, 1, nARR)

    ist = 0
    l   = 1
    r   = nARR

    DO
      IF (r - l < NN) THEN
        DO j = l + 1, r
          tmpIDX = idx1D(j)
          a = arr1D(tmpIDX)
          DO i = j - 1, l, -1
            IF (arr1D(idx1D(i)) <= a) EXIT
            idx1D(i + 1) = idx1D(i)
          END DO
          idx1D(i + 1) = tmpIDX
        END DO

        IF (ist == 0) THEN
          CALL UnsetMessageSource()

          RETURN
        END IF

        r   = stack(ist)
        l   = stack(ist - 1)
        ist = ist - 2
      ELSE
        k = (l + r) / 2

        CALL Swap(idx1D(k), idx1D(l + 1))
        CALL IcompXchg(idx1D(l), idx1D(r))
        CALL IcompXchg(idx1D(l + 1), idx1D(r))
        CALL IcompXchg(idx1D(l), idx1D(l + 1))

        i = l + 1
        j = r
        tmpIDX = idx1D(l + 1)
        a = arr1D(tmpIDX)

        DO
          DO
            i = i + 1
            IF (arr1D(idx1D(i)) > a) EXIT
          END DO

          DO
            j = j - 1
            IF (arr1D(idx1D(j)) < a) EXIT
          END DO

          IF (j < i) EXIT
          CALL Swap(idx1D(i), idx1D(j))
        END DO

        idx1D(l + 1) = idx1D(j)
        idx1D(j) = tmpIDX
        ist = ist + 2

        IF (ist > NSTACK) THEN
          WRITE(tmpStr1, '(a, i0)') 'NSTACK = ', NSTACK
          WRITE(scratchMessage, '(a)') 'The value of the NSTACK parameter is too small: ' // &
                                       TRIM(ADJUSTL(tmpStr1))

          CALL LogMessage(ERROR, scratchMessage)
          CALL UnsetMessageSource()

          IF (PRESENT(status)) status = 2

          RETURN

        END IF

        IF (r - i + 1 >= j - l) THEN
          stack(ist) = r
          stack(ist - 1) = i
          r = j - 1
        ELSE
          stack(ist) = j - 1
          stack(ist - 1) = l
          l = i
        END IF
      END IF
    END DO

    CALL UnsetMessageSource()


    CONTAINS

    SUBROUTINE IcompXchg(i, j)

      IMPLICIT NONE

      ! Global variables
      INTEGER, INTENT(INOUT) :: i, j

      ! Local variables
      INTEGER :: swp

      IF (arr1D(j) < arr1D(i)) THEN
        swp = i
        i   = j
        j   = swp
      END IF

    END SUBROUTINE IcompXchg

  END SUBROUTINE IndexxDouble

  !================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   Q U I C K  S O R T
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Sorts the array arr1D into ascending numerical order using Quicksort.
  !>
  !> @details
  !>   The array arr1D is replaced on output by its sorted rearrangement.
  !>   The parameters NN and NSTACK are defined as:
  !>   - NN is the size of subarrays sorted by straight insertion, and
  !>   - NSTACK is the required auxiliary storage
  !>
  !> @param[in,out]
  !>   arr1D      The one-dimensional array to be sorted
  !> @param[out]
  !>   status     The error status, no error: status = 0 (output)
  !>
  !> @note Adopted from Numerical Recipes for Fortran 90
  !>
  !----------------------------------------------------------------
  SUBROUTINE QuickSort(arr1D, status)

    IMPLICIT NONE

    ! Global variables
    REAL(SZ), DIMENSION(:), INTENT(INOUT) :: arr1D
    INTEGER, OPTIONAL, INTENT(OUT)        :: status

    ! Local variables
    INTEGER, PARAMETER                    :: NN = 15, NSTACK = 50
    REAL(SZ)                              :: a
    INTEGER                               :: nARR
    INTEGER                               :: k, i, j, l, r
    INTEGER                               :: ist, stack(NSTACK)
    CHARACTER(LEN=64)                     :: tmpStr1


    CALL SetMessageSource("QuickSort")

    IF (PRESENT(status)) status = 0

    nARR = size(arr1D, 1)

    ist = 0
    l   = 1
    r   = nARR

    DO
      ! Insertion sort when subarray small enough
      IF (r - l < NN) THEN
        DO j = l + 1, r
          a = arr1D(j)
          DO i = j - 1, l, -1
            IF (arr1D(i) <= a) EXIT
            arr1D(i + 1) = arr1D(i)
          END DO
          arr1D(i + 1) = a
        END DO

        IF (ist == 0) THEN
          CALL UnsetMessageSource()

          RETURN
        END IF

        ! Pop stack and begin a new round of partitioning
        r   = stack(ist)
        l   = stack(ist - 1)
        ist = ist - 2

      ! Choose median of left, center, and right elements as partitioning
      ! element a. Also rearrange so that a(l) <= a(l + 1) <= a(r)
      ELSE
        k = (l + r) / 2

        CALL Swap(arr1D(k), arr1D(l + 1))
        CALL Swap(arr1D(l), arr1D(r), arr1D(l) > arr1D(r))
        CALL Swap(arr1D(l + 1), arr1D(r), arr1D(l + 1) > arr1D(r))
        CALL Swap(arr1D(l), arr1D(l + 1), arr1D(l) > arr1D(l + 1))

        ! Initialize pointers for partitioning
        i = l + 1
        j = r
        a = arr1D(l + 1) ! Partitioning element.

        DO ! Here is the meat.
          ! Scan up to find element >= a
          DO
            i = i + 1
            IF (arr1D(i) > a) EXIT
          END DO

          ! Scan down to find element <= a
          DO
            j = j - 1
            IF (arr1D(j) < a) EXIT
          END DO

          ! Pointers crossed. Exit with partitioning complete.
          IF (j < i) EXIT

          CALL Swap(arr1D(i), arr1D(j)) !Exchange elements.
        END DO

        ! Insert partitioning element
        arr1D(l + 1) = arr1D(j)
        arr1D(j) = a
        ist = ist + 2

        ! Push pointers to larger subarray on stack; process smaller subarray immediately.
        IF (ist > NSTACK) THEN
          WRITE(tmpStr1, '(a, i0)') 'NSTACK = ', NSTACK
          WRITE(scratchMessage, '(a)') 'The value of the NSTACK parameter is too small: ' // &
                                       TRIM(ADJUSTL(tmpStr1))

          CALL LogMessage(ERROR, scratchMessage)
          CALL UnsetMessageSource()

          IF (PRESENT(status)) status = 2

          RETURN

        END IF

        IF (r - i + 1 >= j - l) THEN
          stack(ist)     = r
          stack(ist - 1) = i
          r = j - 1
        ELSE
          stack(ist)     = j - 1
          stack(ist - 1) = l
          l = i
        END IF
      END IF
    END DO

    CALL UnsetMessageSource()

  END SUBROUTINE QuickSort

  !================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   S O R T  2
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Sorts two 1D arrays into ascending numerical order using Quicksort.
  !>
  !> @details
  !>   Sorts the array arr1D into ascending order using Quicksort, while making the corresponding
  !>   rearrangement of the same-size array slv1D. The sorting and rearrangement are performed
  !>   by means of the index array.
  !>
  !> @param[in,out]
  !>   arr1D      The first one-dimensional array to be sorted in ascending order
  !> @param[in,out]
  !>   slv1D      The second one-dimensional array to be sorted in ascending order
  !> @param[out]
  !>   status     The error status, no error: status = 0 (output)
  !>
  !> @note Adopted from Numerical Recipes for Fortran 90
  !>
  !----------------------------------------------------------------
  SUBROUTINE Sort2(arr1D, slv1D, status)

    IMPLICIT NONE

    ! Global variables
    REAL(SZ), DIMENSION(:), INTENT(INOUT) :: arr1D, slv1D
    INTEGER, OPTIONAL, INTENT(OUT)        :: status

    ! Local variables
    INTEGER                         :: nARR, nSLV
    INTEGER, DIMENSION(SIZE(arr1D)) :: idx1D
    CHARACTER(LEN=64)               :: tmpStr1, tmpStr2


    CALL SetMessageSource("Sort2")

    nARR = SIZE(arr1D, 1)
    nSLV = SIZE(slv1D, 1)

    IF (nARR /= nSLV) THEN
      WRITE(tmpStr1, '(a, i0)') 'nARR = ', nARR
      WRITE(tmpStr2, '(a, i0)') 'nSLV = ', nSLV
      WRITE(scratchMessage, '(a)') 'The size of the 1D arrays arr1D and slv1D is not the same: ' // &
                                   TRIM(ADJUSTL(tmpStr1)) // ', ' // TRIM(ADJUSTL(tmpStr2))
      
      CALL LogMessage(ERROR, scratchMessage)
      CALL UnsetMessageSource()

      IF (PRESENT(status)) status = 1

      RETURN
    END IF

    ! Make the index array
    CALL Indexx(arr1D, idx1D, status)
 
    ! Sort the array
    arr1D = arr1D(idx1D)

    ! Rearrange slave
    slv1D = slv1D(idx1D)

    CALL UnsetMessageSource()

  END SUBROUTINE Sort2

  !================================================================================


  !----------------------------------------------------------------
  ! S U B R O U T I N E   A R R A Y  C O P Y  I N T
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Copies the 1D source integer array "src" into the 1D destination array "dest".
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>   src    The one-dimensional array to be copied (integer)
  !> @param[out]
  !>   dest   The copied array (output)
  !> @param[out]
  !>   nCP    The number of elements of "src" array that copied (output)
  !> @param[out]
  !>   nNCP   The number of elements of "src" array that failed to be copied (output)
  !>
  !> @note Adopted from Numerical Recipes for Fortran 90
  !>
  !----------------------------------------------------------------
  SUBROUTINE ArrayCopyInt(src, dest, nCP, nNCP)

    IMPLICIT NONE

    ! Global variables
    INTEGER, DIMENSION(:), INTENT(IN)  :: src
    INTEGER, DIMENSION(:), INTENT(OUT) :: dest
    INTEGER, INTENT(OUT)               :: nCP, nNCP

    nCP  = MIN(SIZE(src), SIZE(dest))
    nNCP = SIZE(src) - nCP
    dest(1:nCP) = src(1:nCP)

  END SUBROUTINE ArrayCopyInt

  !================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   A R R A Y  C O P Y  S I N G L E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Copies the 1D source single precision array "src" into the 1D destination array "dest".
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>   src    The one-dimensional array to be copied (single precision)
  !> @param[out]
  !>   dest   The copied array (output)
  !> @param[out]
  !>   nCP    The number of elements of "src" array that copied (output)
  !> @param[out]
  !>   nNCP   The number of elements of "src" array that failed to be copied (output)
  !>
  !> @note Adopted from Numerical Recipes for Fortran 90
  !>
  !----------------------------------------------------------------
  SUBROUTINE ArrayCopySingle(src, dest, nCP, nNCP)

    IMPLICIT NONE

    ! Global variables
    REAL(SP), DIMENSION(:), INTENT(IN)  :: src
    REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
    INTEGER, INTENT(OUT)                :: nCP, nNCP

    nCP  = MIN(SIZE(src), SIZE(dest))
    nNCP = SIZE(src) - nCP
    dest(1:nCP) = src(1:nCP)

  END SUBROUTINE ArrayCopySingle

  !================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   A R R A Y  C O P Y  D O U B L E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Copies the 1D source double precision array "src" into the 1D destination array "dest".
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>   src    The one-dimensional array to be copied (double precision)
  !> @param[out]
  !>   dest   The copied array (output)
  !> @param[out]
  !>   nCP    The number of elements of "src" array that copied (output)
  !> @param[out]
  !>   nNCP   The number of elements of "src" array that failed to be copied (output)
  !>
  !> @note Adopted from Numerical Recipes for Fortran 90
  !>
  !----------------------------------------------------------------
  SUBROUTINE ArrayCopyDouble(src, dest, nCP, nNCP)

    IMPLICIT NONE

    ! Global variables
    REAL(HP), DIMENSION(:), INTENT(IN)  :: src
    REAL(HP), DIMENSION(:), INTENT(OUT) :: dest
    INTEGER, INTENT(OUT)                :: nCP, nNCP

    nCP  = MIN(SIZE(src), SIZE(dest))
    nNCP = SIZE(src) - nCP
    dest(1:nCP) = src(1:nCP)

  END SUBROUTINE ArrayCopyDouble

  !================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   A R R A Y  E Q U A L  I N T
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Compares two one-dimensional integer arrays for equality.
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>   arr1    The first array in the comparison (integer)
  !> @param[in]
  !>   arr2    The second array in the comparison (integer)
  !>
  !> @return
  !>   myValOut: The value of the comparison (logical). TRUE if all the elements
  !>             of arr1 are equal to all elements of arr2, FALSE otherwise.
  !>
  !----------------------------------------------------------------
  LOGICAL FUNCTION ArrayEqualInt(arr1, arr2) RESULT(myValOut)

    IMPLICIT NONE

    ! Global variables
    INTEGER, DIMENSION(:), INTENT(IN) :: arr1, arr2


    IF (SIZE(arr1) /= SIZE(arr2)) THEN
      myValOut = .FALSE.

      RETURN
    END IF

    myValOut = .TRUE.
    IF (ANY(arr1 - arr2 /= 0)) myValOut = .FALSE.

    RETURN

  END FUNCTION ArrayEqualInt

  !================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   A R R A Y  E Q U A L  S I N G L E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Compares two one-dimensional single precision arrays for equality.
  !>
  !> @details
  !>   The equality is determined using a tolerance of: 0.00000001, such that
  !>   the two arrays are considered to be essentially equal on single
  !>   precision calculations.
  !>
  !> @param[in]
  !>   arr1    The first array in the comparison (single precision)
  !> @param[in]
  !>   arr2    The second array in the comparison (single precision)
  !>
  !> @return
  !>   myValOut: The value of the comparison (logical). TRUE if all the elements
  !>             of arr1 are equal to all elements of arr2, FALSE otherwise.
  !>
  !----------------------------------------------------------------
  LOGICAL FUNCTION ArrayEqualSingle(arr1, arr2) RESULT(myValOut)

    IMPLICIT NONE

    ! Global variables
    REAL(SP), DIMENSION(:), INTENT(IN) :: arr1, arr2

    ! Local variables
    INTEGER :: i


    IF (SIZE(arr1) /= SIZE(arr2)) THEN
      myValOut = .FALSE.

      RETURN
    END IF

    myValOut = .TRUE.
    
    DO i = 1, SIZE(arr1, 1)
      IF (CompareReals(arr1(i), arr2(i), 0.00000001_SP) /= 0) THEN
        myValOut = .FALSE.
        
        EXIT
      END IF
    END DO

    RETURN

  END FUNCTION ArrayEqualSingle

  !================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   A R R A Y  E Q U A L  S I N G L E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Compares two one-dimensional double precision arrays for equality.
  !>
  !> @details
  !>   The equality is determined using a tolerance of: 0.00000001, such that
  !>   the two arrays are considered to be essentially equal on double
  !>   precision calculations.
  !>
  !> @param[in]
  !>   arr1    The first array in the comparison (double precision)
  !> @param[in]
  !>   arr2    The second array in the comparison (double precision)
  !>
  !> @return
  !>   myValOut: The value of the comparison (logical). TRUE if all the elements
  !>             of arr1 are equal to all elements of arr2, FALSE otherwise.
  !>
  !----------------------------------------------------------------
  LOGICAL FUNCTION ArrayEqualDouble(arr1, arr2) RESULT(myValOut)

    IMPLICIT NONE

    ! Global variables
    REAL(HP), DIMENSION(:), INTENT(IN) :: arr1, arr2

    ! Local variables
    INTEGER :: i


    IF (SIZE(arr1) /= SIZE(arr2)) THEN
      myValOut = .FALSE.

      RETURN
    END IF

    myValOut = .TRUE.
    
    DO i = 1, SIZE(arr1, 1)
      IF (CompareReals(arr1(i), arr2(i), 0.0000000000001_HP) /= 0) THEN
        myValOut = .FALSE.
        
        EXIT
      END IF
    END DO

    RETURN

  END FUNCTION ArrayEqualDouble

  !================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   S T R I N G  L E X  C O M P
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Performs a lexical comparison between two strings.
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>   str1          The first string in the comparison
  !> @param[in]
  !>   str2          The second string in the comparison
  !> @param[in]
  !>   mSensitive    Logical flag (.TRUE., .FALSE.) to perform case sensitive lexical comparison
  !>
  !> @return
  !>   myValOut:   The value of the lexical comparison of the two strings (integer)
  !> @verbatim
  !>   myValOut =  0; str1 == str2
  !>   myValOut = -1; str1 < str2
  !>   myValOut =  1; str1 > str2
  !> @endverbatim
  !>
  !----------------------------------------------------------------
  INTEGER FUNCTION StringLexComp(str1, str2, mSensitive) RESULT(myValOut)

    USE PaHM_Utilities, ONLY : ToUpperCase

    IMPLICIT NONE

    ! Global variables
    CHARACTER(LEN=*), INTENT(IN)  :: str1, str2
    LOGICAL, OPTIONAL, INTENT(IN) :: mSensitive

    ! Local variables
    LOGICAL :: sFlag

    sFlag = .TRUE.
    IF (PRESENT(mSensitive)) sFlag = mSensitive

    IF (sFlag) THEN
      IF (TRIM(str1) == TRIM(str2)) THEN
        myValOut = 0
      ELSE IF (TRIM(str1) < TRIM(str2)) THEN
        myValOut = -1
      ELSE
        myValOut = 1
      END IF
    ELSE
      IF (ToUpperCase(TRIM(str1)) == ToUpperCase(TRIM(str2))) THEN
        myValOut = 0
      ELSE IF (ToUpperCase(TRIM(str1)) < ToUpperCase(TRIM(str2))) THEN
        myValOut = -1
      ELSE
        myValOut = 1
      END IF
    END IF

    RETURN

  END FUNCTION StringLexComp

  !================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   S W A P  I N T
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Swaps the contents of a and b (integer).
  !>   increment and a number of terms "n" (including "first").
  !>
  !> @details
  !>   
  !>
  !> @param[in,out]
  !>   a      The first value to be swapped (integer)
  !> @param[in,out]
  !>   b      The second value to be swapped (integer)
  !> @param[in]
  !>   mask   Logical flag to perform the swap, default mask = 'TRUE. (optional)
  !>
  !> @verbatim
  !>   a: The second swapped value
  !>   b: The first swapped value
  !> @endverbatim
  !>
  !> @note Adopted from Numerical Recipes for Fortran 90
  !----------------------------------------------------------------
  SUBROUTINE SwapInt(a, b, mask)

    IMPLICIT NONE

    ! Global variables
    INTEGER, INTENT(INOUT)        :: a, b
    LOGICAL, OPTIONAL, INTENT(IN) :: mask

    ! Local variables
    INTEGER :: dum
    LOGICAL :: mFlag


    mFlag = .TRUE.
    IF (PRESENT(mask)) mFlag = mask

    IF (mFlag) THEN
      dum = a
      a   = b
      b   = dum
    END IF

  END SUBROUTINE SwapInt

  !================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   S W A P  S I N G L E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Swaps the contents of a and b (single precision).
  !>   increment and a number of terms "n" (including "first").
  !>
  !> @details
  !>   
  !>
  !> @param[in,out]
  !>   a      The first value to be swapped (single precision)
  !> @param[in,out]
  !>   b      The second value to be swapped (single precision)
  !> @param[in]
  !>   mask   Logical flag to perform the swap, default mask = 'TRUE. (optional)
  !>
  !> @verbatim
  !>   a: The second swapped value
  !>   b: The first swapped value
  !> @endverbatim
  !>
  !> @note Adopted from Numerical Recipes for Fortran 90
  !----------------------------------------------------------------
  SUBROUTINE SwapSingle(a, b, mask)

    IMPLICIT NONE

    ! Global variables
    REAL(SP), INTENT(INOUT)        :: a, b
    LOGICAL, OPTIONAL, INTENT(IN)  :: mask

    ! Local variables
    REAL(SP) :: dum
    LOGICAL  :: mFlag


    mFlag = .TRUE.
    IF (PRESENT(mask)) mFlag = mask

    IF (mFlag) THEN
      dum = a
      a   = b
      b   = dum
    END IF

  END SUBROUTINE SwapSingle

  !================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   S W A P  D O U B L E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Swaps the contents of a and b (double precision).
  !>   increment and a number of terms "n" (including "first").
  !>
  !> @details
  !>   
  !>
  !> @param[in,out]
  !>   a      The first value to be swapped (double precision)
  !> @param[in]
  !>   b      The second value to be swapped (double precision)
  !> @param[in,out]
  !>   mask   Logical flag to perform the swap, default mask = 'TRUE. (optional)
  !>
  !> @verbatim
  !>   a: The second swapped value
  !>   b: The first swapped value
  !> @endverbatim
  !>
  !> @note Adopted from Numerical Recipes for Fortran 90
  !----------------------------------------------------------------
  SUBROUTINE SwapDouble(a, b, mask)

    IMPLICIT NONE

    ! Global variables
    REAL(HP), INTENT(INOUT)        :: a, b
    LOGICAL, OPTIONAL, INTENT(IN)  :: mask

    ! Local variables
    REAL(HP) :: dum
    LOGICAL  :: mFlag


    mFlag = .TRUE.
    IF (PRESENT(mask)) mFlag = mask

    IF (mFlag) THEN
      dum = a
      a   = b
      b   = dum
    END IF

  END SUBROUTINE SwapDouble

  !================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   S W A P  I N T  V E C
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Swaps the contents of a and b (integer).
  !>   increment and a number of terms "n" (including "first").
  !>
  !> @details
  !>   
  !>
  !> @param[in,out]
  !>   a      The first 1D array to be swapped (integer)
  !> @param[in,out]
  !>   b      The second 1D array to be swapped (integer)
  !> @param[in]
  !>   mask   Logical flag to perform the swap, default mask = 'TRUE. (optional)
  !>
  !> @verbatim
  !>   a: The second swapped 1D array
  !>   b: The first swapped 1D array
  !> @endverbatim
  !>
  !> @note Adopted from Numerical Recipes for Fortran 90
  !----------------------------------------------------------------
  SUBROUTINE SwapIntVec(a, b, mask)

    IMPLICIT NONE

    ! Global variables
    INTEGER, DIMENSION(:), INTENT(INOUT) :: a, b
    LOGICAL, OPTIONAL, INTENT(IN)        :: mask

    ! Local variables
    INTEGER, DIMENSION(SIZE(a)) :: dum
    LOGICAL                     :: mFlag


    mFlag = .TRUE.
    IF (PRESENT(mask)) mFlag = mask

    IF (mFlag) THEN
      dum = a
      a   = b
      b   = dum
    END IF

  END SUBROUTINE SwapIntVec

  !================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   S W A P  S I N G L E  V E C
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Swaps the contents of a and b (single precision).
  !>   increment and a number of terms "n" (including "first").
  !>
  !> @details
  !>   
  !>
  !> @param[in,out]
  !>   a      The first 1D array to be swapped (single precision)
  !> @param[in,out]
  !>   b      The second 1D array to be swapped (single precision)
  !> @param[in]
  !>   mask   Logical flag to perform the swap, default mask = 'TRUE. (optional)
  !>
  !> @verbatim
  !>   a: The second swapped 1D array
  !>   b: The first swapped 1D array
  !> @endverbatim
  !>
  !> @note Adopted from Numerical Recipes for Fortran 90
  !----------------------------------------------------------------
  SUBROUTINE SwapSingleVec(a, b, mask)

    IMPLICIT NONE

    ! Global variables
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a, b
    LOGICAL, OPTIONAL, INTENT(IN)         :: mask

    ! Local variables
    REAL(SP), DIMENSION(SIZE(a)) :: dum
    LOGICAL                      :: mFlag


    mFlag = .TRUE.
    IF (PRESENT(mask)) mFlag = mask

    IF (mFlag) THEN
      dum = a
      a   = b
      b   = dum
    END IF

  END SUBROUTINE SwapSingleVec

  !================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   S W A P  D O U B L E  V E C
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Swaps the contents of a and b (double precision).
  !>   increment and a number of terms "n" (including "first").
  !>
  !> @details
  !>   
  !>
  !> @param[in,out]
  !>   a      The first 1D array to be swapped (double precision)
  !> @param[in,out]
  !>   b      The second 1D array to be swapped (double precision)
  !> @param[in]
  !>   mask   Logical flag to perform the swap, default mask = 'TRUE. (optional)
  !>
  !> @verbatim
  !>   a: The second swapped 1D array
  !>   b: The first swapped 1D array
  !> @endverbatim
  !>
  !> @note Adopted from Numerical Recipes for Fortran 90
  !----------------------------------------------------------------
  SUBROUTINE SwapDoubleVec(a, b, mask)

    IMPLICIT NONE

    ! Global variables
    REAL(HP), DIMENSION(:), INTENT(INOUT) :: a, b
    LOGICAL, OPTIONAL, INTENT(IN)         :: mask

    ! Local variables
    REAL(HP), DIMENSION(SIZE(a)) :: dum
    LOGICAL                      :: mFlag


    mFlag = .TRUE.
    IF (PRESENT(mask)) mFlag = mask

    IF (mFlag) THEN
      dum = a
      a   = b
      b   = dum
    END IF

  END SUBROUTINE SwapDoubleVec

  !================================================================================
 
  !----------------------------------------------------------------
  ! F U N C T I O N   A R T H  I N T
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Returns an arithmetic progression, given a first term "first", an
  !>   increment and a number of terms "n" (including "first").
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>       first    The value of the first term (integer)
  !> @param[in]
  !>   increment    The value of the increment (integer)
  !> @param[in]
  !>           n    The total number of terms in the return 1D array (integer)
  !>
  !> @return
  !>   arthOut:     The 1D array that contains the arithmetic progression (integer)
  !>
  !> @note Adopted from Numerical Recipes for Fortran 90
  !----------------------------------------------------------------
  pure FUNCTION ArthInt(first, increment, n) RESULT(arthOut)

    IMPLICIT NONE

    ! Global variables
    INTEGER, INTENT(IN)   :: first, increment
    INTEGER, INTENT(IN)   :: n
    INTEGER, DIMENSION(n) :: arthOut

    ! Local variables
    INTEGER, PARAMETER :: NPARTH = 16, NPARTH2 = 8
    INTEGER :: k, k2
    INTEGER :: temp


    IF (n > 0) arthOut(1) = first

    IF (n <= NPARTH) THEN
      DO k = 2, n
        arthOut(k) = arthOut(k - 1) + increment
      END DO
    ELSE
      DO k = 2, NPARTH2
        arthOut(k) = arthOut(k - 1) + increment
      END DO

      temp = increment * NPARTH2
      k = NPARTH2

      DO
        IF (k >= n) EXIT
        k2 = k + k
        arthOut(k + 1:min(k2, n)) = temp + arthOut(1:min(k, n - k))
        temp = temp + temp
        k = k2
      END DO
    END IF

    RETURN

  END FUNCTION ArthInt

  !================================================================================ 

  !----------------------------------------------------------------
  ! F U N C T I O N   A R T H  S I N G L E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Returns an arithmetic progression, given a first term "first", an
  !>   increment and a number of terms "n" (including "first").
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>       first    The value of the first term (single precision)
  !> @param[in]
  !>   increment    The value of the increment (single precision)
  !> @param[in]
  !>           n    The total number of terms in the return 1D array (integer)
  !>
  !> @return
  !>   arthOut:     The 1D array that contains the arithmetic progression (single precision)
  !>
  !> @note Adopted from Numerical Recipes for Fortran 90
  !----------------------------------------------------------------
  pure FUNCTION ArthSingle(first, increment, n) RESULT(arthOut)

    IMPLICIT NONE

    ! Global variables
    REAL(SP), INTENT(IN)   :: first, increment
    INTEGER, INTENT(IN)    :: n
    REAL(SP), DIMENSION(n) :: arthOut

    ! Local variables
    INTEGER, PARAMETER :: NPARTH = 16, NPARTH2 = 8
    INTEGER  :: k, k2
    REAL(SP) :: temp


    IF (n > 0) arthOut(1) = first

    IF (n <= NPARTH) THEN
      DO k = 2, n
        arthOut(k) = arthOut(k - 1) + increment
      END DO
    ELSE
      DO k = 2, NPARTH2
        arthOut(k) = arthOut(k - 1) + increment
      END DO

      temp = increment * NPARTH2
      k = NPARTH2

      DO
        IF (k >= n) EXIT
        k2 = k + k
        arthOut(k + 1:min(k2, n)) = temp + arthOut(1:min(k, n - k))
        temp = temp + temp
        k = k2
      END DO
    END IF

    RETURN

  END FUNCTION ArthSingle

  !================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   A R T H  D O U B L E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Returns an arithmetic progression, given a first term "first", an
  !>   increment and a number of terms "n" (including "first").
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>       first    The value of the first term (double precision)
  !> @param[in]
  !>   increment    The value of the increment (double precision)
  !> @param[in]
  !>           n    The total number of terms in the return 1D array (integer)
  !>
  !> @return
  !>   arthOut:     The 1D array that contains the arithmetic progression (double precision)
  !>
  !> @note Adopted from Numerical Recipes for Fortran 90
  !----------------------------------------------------------------
  pure FUNCTION ArthDouble(first, increment, n) RESULT(arthOut)

    IMPLICIT NONE

    ! Global variables
    REAL(HP), INTENT(IN)   :: first, increment
    INTEGER, INTENT(IN)    :: n
    REAL(HP), DIMENSION(n) :: arthOut

    ! Local variables
    INTEGER, PARAMETER :: NPARTH = 16, NPARTH2 = 8
    INTEGER  :: k, k2
    REAL(HP) :: temp


    IF (n > 0) arthOut(1) = first

    IF (n <= NPARTH) THEN
      DO k = 2, n
        arthOut(k) = arthOut(k - 1) + increment
      END DO
    ELSE
      DO k = 2, NPARTH2
        arthOut(k) = arthOut(k - 1) + increment
      END DO

      temp = increment * NPARTH2
      k = NPARTH2

      DO
        IF (k >= n) EXIT
        k2 = k + k
        arthOut(k + 1:min(k2, n)) = temp + arthOut(1:min(k, n - k))
        temp = temp + temp
        k = k2
      END DO
    END IF

    RETURN

  END FUNCTION ArthDouble

  !================================================================================ 

END MODULE SortUtils
