      REAL FUNCTION AKI (OM, BETA)

! ----------------------------------------------------------------------

!**** *AKI* - FUNCTION TO COMPUTE WAVE NUMBER.

!     G. KOMEN, P. JANSSEN   KNMI        01/06/1986

!*    PURPOSE.
!     -------

!       *AKI* COMPUTES THE WAVE NUMBER AS FUNCTION OF
!             CIRCULAR FREQUENCY AND WATER DEPTH.

!**   INTERFACE.
!     ----------

!       *FUNCTION* *AKI (OM, BETA)*
!          *OM*      - CIRCULAR FREQUENCY.
!          *BETA*    - WATER DEPTH.

!     METHOD.
!     -------

!       NEWTONS METHOD TO SOLVE THE DISPERSION RELATION IN SHALLOW
!       WATER.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------


!      USE YOWPCONS , ONLY : G     ,DKMAX
       USE DATAPOOL, ONLY : G => G9, DKMAX 

!*    *PARAMETER*  RELATIVE ERROR LIMIT OF NEWTON'S METHOD.

      REAL, PARAMETER :: EBS = 0.0001

! ----------------------------------------------------------------------

!*    1. START VALUE:  MAXIMUM FROM DEEP  AND EXTREM SHALLOW WATER
!                      WAVE NUMBER.
!        ---------------------------------------------------------

      AKM1=OM**2/(4.*G)
      AKM2=OM/(2.*SQRT(G*BETA))
      AO=MAX(AKM1,AKM2)

! ----------------------------------------------------------------------

!*    2. ITERATION LOOP.
!        ---------------

 2000 CONTINUE
      AKP = AO
      BO = BETA*AO
      IF (BO.GT.DKMAX) THEN
        AKI = OM**2/G
      ELSE
        TH = G*AO*TANH(BO)
        STH = SQRT(TH)
        AO = AO+(OM-STH)*STH*2./(TH/AO+G*BO/COSH(BO)**2)
        IF (ABS(AKP-AO).GT.EBS*AO) GO TO 2000
        AKI = AO
      ENDIF

      RETURN
      END




