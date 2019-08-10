#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_SNL4_WW3
      USE DATAPOOL
!
!
! this code is based on Hendrik Tolmans code from WWIII v.1.18, 1999
!
!
      IMPLICIT NONE
      INTEGER                 :: IFR, ITH, ISP, ITHP, ITHP1, ITHM,ITHM1,IFRP, IFRP1, IFRM, IFRM1

      INTEGER, ALLOCATABLE    :: IF1(:), IF2(:), IF3(:), IF4(:),IF5(:), IF6(:), IF7(:), IF8(:)
                                 IT1(:), IT2(:), IT3(:), IT4(:),IT5(:), IT6(:), IT7(:), IT8(:)

      REAL                    :: DELTH3, DELTH4, LAMM2, LAMP2, CTHP, WTHP, WTHP1, CTHM, WTHM, WTHM1
      REAL(rkind)             :: XFRLN, WFRP, WFRP1, WFRM, WFRM1, FR, AF11A
!
      NFR     = NK
!
! 1.  Internal angles of quadruplet.
!
      LAMM2  = (1.-LAM)**2
      LAMP2  = (1.+LAM)**2
      DELTH3 = ACOS( (LAMM2**2+4.-LAMP2**2) / (4.*LAMM2) )
      DELTH4 = ASIN(-SIN(DELTH3)*LAMM2/LAMP2)
!
! 2.  Lambda dependend weight factors.
!
      DAL1   = 1. / (1.+LAM)**4
      DAL2   = 1. / (1.-LAM)**4
      DAL3   = 2. * DAL1 * DAL2
!
! 3.  Directional indices.
!
      CTHP   = ABS(DELTH4/DTH)
      ITHP   = INT(CTHP)
      ITHP1  = ITHP + 1
      WTHP   = CTHP - REAL(ITHP)
      WTHP1  = 1.- WTHP
!
      CTHM   = ABS(DELTH3/DTH)
      ITHM   = INT(CTHM)
      ITHM1  = ITHM + 1
      WTHM   = CTHM - REAL(ITHM)
      WTHM1  = 1.- WTHM
!
! 4.  Frequency indices.
!
      XFRLN  = LOG(XFR)
!
      IFRP   = INT( LOG(1.+LAM) / XFRLN )
      IFRP1  = IFRP + 1
      WFRP   = (1.+LAM - XFR**IFRP) / (XFR**IFRP1 - XFR**IFRP)
      WFRP1  = 1. - WFRP

      IFRM   = INT( LOG(1.-LAM) / XFRLN )
      IFRM1  = IFRM - 1
      WFRM   = (XFR**IFRM -(1.-LAM)) / (XFR**IFRM - XFR**IFRM1)
      WFRM1  = 1. - WFRM
!
! 5.  Range of calculations
!
      NFRHGH = NFR + IFRP1 - IFRM1
      NFRCHG = NFR - IFRM1
      NSPECY = NFRHGH * NTH
      NSPECX = NFRCHG * NTH
!
! 6.  Allocate arrays or check array sizes
!
!
! 6.  Allocate arrays or check array sizes
!
      ALLOCATE ( IP11(NSPECX),IP12(NSPECX),IP13(NSPECX),IP14(NSPECX), &
                 IM11(NSPECX),IM12(NSPECX),IM13(NSPECX),IM14(NSPECX), &
                 IP21(NSPECX),IP22(NSPECX),IP23(NSPECX),IP24(NSPECX), &
                 IM21(NSPECX),IM22(NSPECX),IM23(NSPECX),IM24(NSPECX), &
                 IC11(NSPEC), IC12(NSPEC), IC21(NSPEC), IC22(NSPEC),  &
                 IC31(NSPEC), IC32(NSPEC), IC41(NSPEC), IC42(NSPEC),  &
                 IC51(NSPEC), IC52(NSPEC), IC61(NSPEC), IC62(NSPEC),  &
                 IC71(NSPEC), IC72(NSPEC), IC81(NSPEC), IC82(NSPEC),  &
                 AF11(NSPECX) )
!
      ALLOCATE ( IF1(NFRCHG), IF2(NFRCHG), IF3(NFRCHG), IF4(NFRCHG),  &
                 IF5(NFRCHG), IF6(NFRCHG), IF7(NFRCHG), IF8(NFRCHG),  &
                 IT1(NTH), IT2(NTH), IT3(NTH), IT4(NTH),              &
                 IT5(NTH), IT6(NTH), IT7(NTH), IT8(NTH) )
!
! 7.  Spectral addresses
!
      DO IFR=1, NFRCHG
        IF1(IFR) =           IFR+IFRP
        IF2(IFR) =           IFR+IFRP1
        IF3(IFR) = MAX ( 0 , IFR+IFRM  )
        IF4(IFR) = MAX ( 0 , IFR+IFRM1 )
        IF5(IFR) = MAX ( 0 , IFR-IFRP  )
        IF6(IFR) = MAX ( 0 , IFR-IFRP1 )
        IF7(IFR) =           IFR-IFRM
        IF8(IFR) =           IFR-IFRM1
        END DO
!
      DO ITH=1, NTH
        IT1(ITH) = ITH + ITHP
        IT2(ITH) = ITH + ITHP1
        IT3(ITH) = ITH + ITHM
        IT4(ITH) = ITH + ITHM1
        IT5(ITH) = ITH - ITHP
        IT6(ITH) = ITH - ITHP1
        IT7(ITH) = ITH - ITHM
        IT8(ITH) = ITH - ITHM1
        IF ( IT1(ITH).GT.NTH) IT1(ITH) = IT1(ITH) - NTH
        IF ( IT2(ITH).GT.NTH) IT2(ITH) = IT2(ITH) - NTH
        IF ( IT3(ITH).GT.NTH) IT3(ITH) = IT3(ITH) - NTH
        IF ( IT4(ITH).GT.NTH) IT4(ITH) = IT4(ITH) - NTH
        IF ( IT5(ITH).LT. 1 ) IT5(ITH) = IT5(ITH) + NTH
        IF ( IT6(ITH).LT. 1 ) IT6(ITH) = IT6(ITH) + NTH
        IF ( IT7(ITH).LT. 1 ) IT7(ITH) = IT7(ITH) + NTH
        IF ( IT8(ITH).LT. 1 ) IT8(ITH) = IT8(ITH) + NTH
        END DO
!
      DO ISP=1, NSPECX
        IFR       = 1 + (ISP-1)/NTH
        ITH       = 1 + MOD(ISP-1,NTH)
        IP11(ISP) = IT2(ITH) + (IF2(IFR)-1)*NTH
        IP12(ISP) = IT1(ITH) + (IF2(IFR)-1)*NTH
        IP13(ISP) = IT2(ITH) + (IF1(IFR)-1)*NTH
        IP14(ISP) = IT1(ITH) + (IF1(IFR)-1)*NTH
        IM11(ISP) = IT8(ITH) + (IF4(IFR)-1)*NTH
        IM12(ISP) = IT7(ITH) + (IF4(IFR)-1)*NTH
        IM13(ISP) = IT8(ITH) + (IF3(IFR)-1)*NTH
        IM14(ISP) = IT7(ITH) + (IF3(IFR)-1)*NTH
        IP21(ISP) = IT6(ITH) + (IF2(IFR)-1)*NTH
        IP22(ISP) = IT5(ITH) + (IF2(IFR)-1)*NTH
        IP23(ISP) = IT6(ITH) + (IF1(IFR)-1)*NTH
        IP24(ISP) = IT5(ITH) + (IF1(IFR)-1)*NTH
        IM21(ISP) = IT4(ITH) + (IF4(IFR)-1)*NTH
        IM22(ISP) = IT3(ITH) + (IF4(IFR)-1)*NTH
        IM23(ISP) = IT4(ITH) + (IF3(IFR)-1)*NTH
        IM24(ISP) = IT3(ITH) + (IF3(IFR)-1)*NTH
        END DO
      DO ISP=1, NSPEC
        IFR       = 1 + (ISP-1)/NTH
        ITH       = 1 + MOD(ISP-1,NTH)
        IC11(ISP) = IT6(ITH) + (IF6(IFR)-1)*NTH
        IC21(ISP) = IT5(ITH) + (IF6(IFR)-1)*NTH
        IC31(ISP) = IT6(ITH) + (IF5(IFR)-1)*NTH
        IC41(ISP) = IT5(ITH) + (IF5(IFR)-1)*NTH
        IC51(ISP) = IT4(ITH) + (IF8(IFR)-1)*NTH
        IC61(ISP) = IT3(ITH) + (IF8(IFR)-1)*NTH
        IC71(ISP) = IT4(ITH) + (IF7(IFR)-1)*NTH
        IC81(ISP) = IT3(ITH) + (IF7(IFR)-1)*NTH
        IC12(ISP) = IT2(ITH) + (IF6(IFR)-1)*NTH
        IC22(ISP) = IT1(ITH) + (IF6(IFR)-1)*NTH
        IC32(ISP) = IT2(ITH) + (IF5(IFR)-1)*NTH
        IC42(ISP) = IT1(ITH) + (IF5(IFR)-1)*NTH
        IC52(ISP) = IT8(ITH) + (IF8(IFR)-1)*NTH
        IC62(ISP) = IT7(ITH) + (IF8(IFR)-1)*NTH
        IC72(ISP) = IT8(ITH) + (IF7(IFR)-1)*NTH
        IC82(ISP) = IT7(ITH) + (IF7(IFR)-1)*NTH
        END DO
!
      DEALLOCATE ( IF1, IF2, IF3, IF4, IF5, IF6, IF7, IF8,  &
                   IT1, IT2, IT3, IT4, IT5, IT6, IT7, IT8 )
!
! 8.  Fill scaling array (f**11)
!
      DO IFR=1, NFR
        AF11A  = (SIG(IFR)*TPIINV)**11
        DO ITH=1, NTH
          AF11(ITH+(IFR-1)*NTH) = AF11A
          END DO
        END DO
!
      FR     = SIG(NFR)*TPIINV
      DO IFR=NFR+1, NFRCHG
        FR     = FR * XFR
        AF11A  = FR**11
        DO ITH=1, NTH
          AF11(ITH+(IFR-1)*NTH) = AF11A
          END DO
        END DO
!
! 9.  Interpolation weights
!
      AWG1   = WTHP  * WFRP
      AWG2   = WTHP1 * WFRP
      AWG3   = WTHP  * WFRP1
      AWG4   = WTHP1 * WFRP1
      AWG5   = WTHM  * WFRM
      AWG6   = WTHM1 * WFRM
      AWG7   = WTHM  * WFRM1
      AWG8   = WTHM1 * WFRM1
!
      SWG1   = AWG1**2
      SWG2   = AWG2**2
      SWG3   = AWG3**2
      SWG4   = AWG4**2
      SWG5   = AWG5**2
      SWG6   = AWG6**2
      SWG7   = AWG7**2
      SWG8   = AWG8**2
!/
      END SUBROUTINE INSNL1
!**********************************************************************
!*                                                                    *
!**********************************************************************

