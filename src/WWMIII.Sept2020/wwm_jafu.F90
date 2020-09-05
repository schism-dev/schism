      INTEGER FUNCTION JAFU (CL, J, IAN)
      USE DATAPOOL, ONLY : RKIND
      IMPLICIT NONE

      REAL(rkind), INTENT(IN) :: CL
      INTEGER, INTENT(IN) :: J, IAN
 
      INTEGER             :: IDPH, JA 

! ----------------------------------------------------------------------

!**** *JAFU* - FUNCTION TO COMPUTE THE INDEX ARRAY FOR THE
!              ANGLES OF THE INTERACTING WAVENUMBERS.

!     S. HASSELMANN        MPIFM        01/12/1985.

!*    PURPOSE.
!     --------

!       INDICES DEFINING BINS IN FREQUENCY AND DIRECTION PLANE INTO
!       WHICH NONLINEAR ENERGY TRANSFER INCREMENTS ARE STORED. NEEDED
!       FOR COMPUTATION OF THE NONLINEAR ENERGY TRANSFER.

!**   INTERFACE.
!     ----------

!       *FUNCTION* *JAFU (CL, J, IAN)*
!          *CL*  - WEIGHTS.
!          *J*   - INDEX IN ANGULAR ARRAY.
!          *IAN* - NUMBER OF ANGLES IN ARRAY.

!     METHOD.
!     -------

!       SEE REFERENCE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!        S. HASSELMANN AND K. HASSELMANN,JPO, 1985 B.

! ----------------------------------------------------------------------

      IDPH = CL
      JA = J+IDPH
      IF (JA.LE.0)   JA = IAN+JA-1
      IF (JA.GE.IAN) JA = JA-IAN+1
      JAFU = JA

      RETURN
      END




