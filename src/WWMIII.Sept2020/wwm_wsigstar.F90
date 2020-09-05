      SUBROUTINE WSIGSTAR (IJS, IJL, USNEW, Z0NEW, WSTAR, SIG_N)
! ----------------------------------------------------------------------

!**** *WSIGSTAR* - COMPUTATION OF STANDARD DEVIATION OF USTAR.

!*    PURPOSE.
!     ---------

!     COMPUTES THE STANDARD DEVIATION OF USTAR DUE TO SMALL SCALE GUSTINESS.

!**   INTERFACE.
!     ----------

!     *CALL* *WSIGSTAR (IJS, IJL, USNEW, Z0NEW, WSTAR, SIG_N)
!             *IJS*   - INDEX OF FIRST GRIDPOINT.
!             *IJL*   - INDEX OF LAST GRIDPOINT.
!             *USNEW* - NEW FRICTION VELOCITY IN M/S.
!             *Z0NEW* - ROUGHNESS LENGTH IN M.
!             *WSTAR* - FREE CONVECTION VELOCITY SCALE (M/S).
!             *SIG_N* - ESTINATED STANDARD DEVIATION OF USTAR.


!     METHOD.
!     -------

!     USE PANOFSKY (1991) TO EXPRESS THE STANDARD DEVIATION OF U10 IN TERMS
!     USTAR AND (Zi/L) THE MONIN-OBOKHOV LENGTH (Zi THE INVERSION HEIGHT).
!     (but with the background gustiness set to 0.)
!     and USTAR=SQRT(Cd)*U10 to DERIVE THE STANDARD DEVIATION OF USTAR.
!     WITH CD=A+B*U10 (see below).

!     EXTERNALS.
!     ----------

!       NONE.

!     MODIFICATIONS
!     -------------

!     REFERENCE.
!     ----------

!     SEE SECTION 3.2.1 OF THE WAM DOCUMENTATION.

! ----------------------------------------------------------------------

!      USE YOWCOUP  , ONLY : XKAPPA
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

      USE DATAPOOL, ONLY : XKAPPA, RKIND

! ----------------------------------------------------------------------
      IMPLICIT NONE

      REAL, PARAMETER :: A = 0.8/1000.
      REAL, PARAMETER :: B = 0.08/1000.

      INTEGER :: IJ,IJS,IJL

      REAL(rkind) :: ONETHIRD, BG_GUST
      REAL(rkind) :: U10, C_D, DC_DDU, SIG_CONV
      REAL(rkind) :: ZHOOK_HANDLE
      REAL(rkind), DIMENSION(IJS:IJL) :: USNEW, Z0NEW, WSTAR
      REAL(rkind), DIMENSION(IJS:IJL) :: SIG_N

      !IF (LHOOK) CALL DR_HOOK('WSIGSTAR',0,ZHOOK_HANDLE)

! ----------------------------------------------------------------------

      ONETHIRD = 1./3.
      BG_GUST  = 0.        ! NO BACKGROUND GUSTINESS (S0 12. IS NOT USED)

      DO IJ=IJS,IJL
!
!       IN THE FOLLOWING U10 IS ESTIMATED ASSUMING EVERYTHING IS
!       BASED ON U*
!
        U10 = USNEW(IJ)/XKAPPA*LOG(10./Z0NEW(IJ))
        C_D = A+B*U10
        DC_DDU = B
        SIG_CONV = 1. + 0.5*U10/C_D*DC_DDU
        SIG_N(IJ) = MIN(0.5, SIG_CONV * &
     &                        (BG_GUST*USNEW(IJ)**3+ &
     &                         0.5*XKAPPA*WSTAR(IJ)**3)**ONETHIRD &
     &                           /U10 &
     &                 )
      ENDDO

      !IF (LHOOK) CALL DR_HOOK('WSIGSTAR',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE WSIGSTAR
