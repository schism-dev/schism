#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GRADDEP()
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER     :: IP
         REAL(rkind) :: GDL, GDD

         DDEP(:,:) = 0.0

         SELECT CASE (DIMMODE)

            CASE (1)
#ifdef MPI_PARALL_GRID
               call parallel_abort('WWM - 1d mode cannot work with SCHISM ')
#endif
               CALL DIFFERENTIATE_XDIR(DEP,DDEP(:,1))

               IF (LSPHE) THEN
                  DO IP = 1, MNP
                     DDEP(IP,1) = DDEP(IP,1)/( DEGRAD*REARTH*COS(YP(IP)*DEGRAD) )
                  END DO
               END IF

               IF (LTEST .AND. ITEST > 100) THEN
                  WRITE(STAT%FHNDL,*) 'Gradients of depth '
                  WRITE(STAT%FHNDL,*) '@D/@X '
                  DO IP = 1, MNP
                     WRITE(STAT%FHNDL,'(1X,I5,3F10.5)') IP, DDEP(IP,1)
                  END DO
                  FLUSH(STAT%FHNDL)
               END IF

            CASE (2)

               CALL DIFFERENTIATE_XYDIR(DEP,DDEP(:,1),DDEP(:,2))

               IF (LSPHE) THEN
                  DO IP = 1, MNP
                     DDEP(IP,1) = DDEP(IP,1)/( DEGRAD*REARTH*COS(YP(IP)*DEGRAD) )
                     DDEP(IP,2) = DDEP(IP,2)/( DEGRAD*REARTH )
                  END DO
               END IF

               IF (LSLOP) THEN
                  DO IP = 1, MNP
                     GDL = SQRT(DDEP(IP,1)**2 + DDEP(IP,2)**2) ! Achtung Gradient mit SQRT berechnet ...
                     GDD = MyATAN2(DDEP(IP,2), DDEP(IP,1))
                     IF (GDL < 1.0E-8) CYCLE
                     GDL = SQRT(GDL)
                     IF (GDL > SLMAX) THEN
                        DDEP(IP,1) = SLMAX*COS(GDD)
                        DDEP(IP,2) = SLMAX*SIN(GDD)
                        WRITE(STAT%FHNDL,*) IP, SLMAX, GDL, GDD , 'MAXSLOPE'
                     END IF
                  END DO
                  FLUSH(STAT%FHNDL)
               END IF

            CASE DEFAULT
         END SELECT

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef TIMINGS
      SUBROUTINE WAV_MY_WTIME(wtime)
      USE DATAPOOL, only : rkind
      implicit none
      real(rkind), intent(inout) :: wtime
# ifdef MPI_PARALL_GRID
      real(8) mpi_wtime
      wtime=mpi_wtime()
# else
      CALL cpu_time(wtime)
# endif
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GRADCURT()
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP
         DCUX(:,:) = 0.0
         DCUY(:,:) = 0.0

         SELECT CASE (DIMMODE)
            CASE (1)
               CALL DIFFERENTIATE_XDIR(CURTXY(:,1),DCUX(:,1))
               CALL DIFFERENTIATE_XDIR(CURTXY(:,2),DCUY(:,1))
               IF (LSPHE) THEN
                  DO IP = 1, MNP
                     DCUX(IP,1) = DCUX(IP,1)/( DEGRAD*REARTH*COS(YP(IP)*DEGRAD) )
                     DCUY(IP,1) = DCUY(IP,1)/( DEGRAD*REARTH*COS(YP(IP)*DEGRAD) )
                  END DO
               END IF
               IF (LTEST .AND. ITEST > 100) THEN
                  WRITE(STAT%FHNDL,*) 'Gradients of depth and current'
                  WRITE(STAT%FHNDL,*) '@U/@X     @V/@X'
                  DO IP = 1, MNP
                     WRITE(STAT%FHNDL,'(1X,I5,3F10.5)') IP, DCUX, DCUY
                  END DO
                  FLUSH(STAT%FHNDL)
               END IF
            CASE (2)
               CALL DIFFERENTIATE_XYDIR(CURTXY(:,1),DCUX(:,1),DCUX(:,2))
               CALL DIFFERENTIATE_XYDIR(CURTXY(:,2),DCUY(:,1),DCUY(:,2))
               IF (LSPHE) THEN
                  DO IP = 1, MNP
                     DCUX(IP,1) = DCUX(IP,1)/( DEGRAD*REARTH*COS(YP(IP)*DEGRAD) )
                     DCUY(IP,1) = DCUY(IP,1)/( DEGRAD*REARTH*COS(YP(IP)*DEGRAD) )
                     DCUX(IP,2) = DCUX(IP,2)/( DEGRAD*REARTH )
                     DCUY(IP,2) = DCUY(IP,2)/( DEGRAD*REARTH )
                  END DO
               END IF
               IF (LTEST .AND. ITEST > 100) THEN
                  WRITE(STAT%FHNDL,*) 'The Gradient of Depth and Current'
                  WRITE(STAT%FHNDL,*) ' @U/@X    @U/@Y    @V/@X    @V/@Y'
                  DO IP = 1, MNP
                     WRITE(STAT%FHNDL,'(1X,I5,4F15.7)') IP, DCUX(IP,1), DCUX(IP,2), DCUY(IP,1), DCUY(IP,2)
                  END DO
                  FLUSH(STAT%FHNDL)
               END IF
            CASE DEFAULT
         END SELECT

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GRAD_CG_K()
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IS

         DCGDX(:,:) = 0.0
         DCGDY(:,:) = 0.0
         DWKDX(:,:) = 0.0
         DWKDY(:,:) = 0.0

         SELECT CASE (DIMMODE)
            CASE (1)
              DO IS = 1, MSC
                CALL DIFFERENTIATE_XDIR(CG(IS,:),DCGDX(:,IS))
                CALL DIFFERENTIATE_XDIR(WK(IS,:),DWKDX(:,IS))
              ENDDO
            CASE (2)
              DO IS = 1, MSC
                CALL DIFFERENTIATE_XYDIR(CG(IS,:),DCGDX(:,IS),DCGDY(:,IS))
                CALL DIFFERENTIATE_XYDIR(WK(IS,:),DWKDX(:,IS),DWKDY(:,IS))
              ENDDO
            CASE DEFAULT
         END SELECT

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SMOOTH( BETA, MNP, XP, YP, VAR )
         USE DATAPOOL, ONLY : RKIND, ZERO, ONE, TWO
         IMPLICIT NONE

         INTEGER :: MNP
         REAL(rkind), INTENT(IN)    :: XP(MNP), YP(MNP)
         REAL(rkind), INTENT(INOUT) :: VAR(MNP)
         REAL(rkind), INTENT(IN)    :: BETA
         REAL(rkind)                :: VART(MNP)
         REAL(rkind)                :: SW, SWQ, DISX, DISY, DIST, DIS
         INTEGER                    :: I, J
         
         DO I = 1, MNP
            SW = ZERO
            SWQ = ZERO
            DO J = 1, MNP               
               DISX = (XP(I) - XP(J))**2
               DISY = (YP(I) - YP(J))**2
               DIST = DISX + DISY
               IF (DIST > TINY(1.)) THEN
                  DIS = SQRT(DIST)**BETA
               ELSE
                  DIS = ONE
               END IF 
               SW = SW + DIS
               SWQ = SWQ + DIS*VAR(J)
            END DO
            VART(I) = SWQ / SW
         END DO
         VAR(:) = VART(:)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SMOOTH_V2(VAR_IN, VAR_OUT)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(IN)  :: VAR_IN(MNP)
      REAL(rkind), INTENT(OUT) :: VAR_OUT(MNP)
      INTEGER IP, IADJ, IP_ADJ
      REAL(rkind) :: SumVAR, SumSI, eVal
      DO IP = 1, NP_RES
        SumVAR=SI(IP) * VAR_IN(IP)
        SumSI =SI(IP)
        DO IADJ=1,VERT_DEG(IP)
          IP_ADJ=LIST_ADJ_VERT(IADJ,IP)
          SumVAR=SumVAR + SI(IP_ADJ)*VAR_IN(IP_ADJ)
          SumSI =SumSI  + SI(IP_ADJ)
        END DO
        eVal=SumVAR/SumSI
        VAR_OUT(IP)=eVal
      END DO
#ifdef MPI_PARALL_GRID
      CALL exchange_p2d(VAR_OUT)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SMOOTH_ON_TRIANGLE(VAR_IN, VAR_OUT, nbIter)
      USE DATAPOOL
      IMPLICIT NONE
      integer, intent(in) :: nbIter
      REAL(rkind), INTENT(IN)  :: VAR_IN(MNE)
      REAL(rkind), INTENT(OUT) :: VAR_OUT(MNE)
      INTEGER IE, IEadj, I, nb, iIter
      REAL(rkind) :: eVal, eValB
      REAL(rkind) :: VARextent(MNEextent)
      VAR_OUT=VAR_IN
      DO iIter=1,nbIter
        VARextent(1:MNE)=VAR_OUT
#ifdef MPI_PARALL_GRID
        CALL TRIG_SYNCHRONIZATION(VARextent)
#endif
        DO IE=1,MNE
          eVal=ZERO
          nb=0
          DO I=1,3
            IEadj=IEneighbor(I,IE)
            IF (IEadj .gt. 0) THEN
              eVal=eVal + VARextent(IEadj)
              nb=nb+1
            END IF
          END DO
          eValB=(MyREAL(6-nb)*VARextent(IE) + eVal)/6.0_rkind
          VAR_OUT(IE)=eValB
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DIFFERENTIATE_XDIR(VAR, DVDX)
         USE DATAPOOL
         IMPLICIT NONE
         REAL(rkind), INTENT(IN)  :: VAR(MNP)
         REAL(rkind), INTENT(OUT) :: DVDX(MNP)
         INTEGER           :: IP
         REAL(rkind)       :: TMP1, TMP2

         DVDX(1)   = (VAR(2)-VAR(1))/(XP(2)-XP(1))
         DVDX(MNP) = (VAR(MNP)-VAR(MNP-1))/(XP(MNP)-XP(MNP-1))
         DO IP = 2, MNP-1
!             DVDX(IP) = ( VAR(IP)-VAR(IP-1) ) / ( XP(IP)-XP(IP-1) )
!             DVDX(IP) = (VAR(IP+1)-VAR(IP))/(XP(IP+1)-XP(IP))
!             TMP1     = (VAR(IP)-VAR(IP-1))/(XP(IP)-XP(IP-1))
!             TMP2     = (VAR(IP+1)-VAR(IP))/(XP(IP+1)-XP(IP))
!             DVDX(IP) = 0.5 * (TMP1 + TMP2)
             TMP1     = VAR(IP+1)-VAR(IP-1)
             TMP2     = XP(IP+1)-XP(IP-1)
             DVDX(IP) = TMP1/TMP2
         END DO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DIFFERENTIATE_XYDIR(VAR, DVDX, DVDY)
         USE DATAPOOL
         IMPLICIT NONE
         REAL(rkind), INTENT(IN)  :: VAR(MNP)
         REAL(rkind), INTENT(OUT) :: DVDX(MNP), DVDY(MNP)
         INTEGER           :: NI(3)
         INTEGER           :: IE, I1, I2, I3, IP
         REAL(rkind)            :: DEDY(3),DEDX(3)
         REAL(rkind)            :: DVDXIE, DVDYIE

         REAL(rkind)            :: WEI(MNP)

         WEI(:)  = 0.0_rkind
         DVDX(:) = 0.0_rkind
         DVDY(:) = 0.0_rkind

         DO IE = 1, MNE 
            NI = INE(:,IE)
            I1 = INE(1,IE)
            I2 = INE(2,IE)
            I3 = INE(3,IE)
            WEI(NI) = WEI(NI) + 2.*TRIA(IE)
            DEDX(1) = IEN(1,IE)
            DEDX(2) = IEN(3,IE)
            DEDX(3) = IEN(5,IE)
            DEDY(1) = IEN(2,IE)
            DEDY(2) = IEN(4,IE)
            DEDY(3) = IEN(6,IE)
            DVDXIE  = DOT_PRODUCT( VAR(NI),DEDX)
            DVDYIE  = DOT_PRODUCT( VAR(NI),DEDY)
            DVDX(NI) = DVDX(NI) + DVDXIE
            DVDY(NI) = DVDY(NI) + DVDYIE
         END DO

         DVDX(:) = DVDX(:)/WEI(:)
         DVDY(:) = DVDY(:)/WEI(:)

         DO IP = 1, MNP
           IF (DEP(IP) .LT. DMIN) THEN
             DVDX(IP) = 0.
             DVDY(IP) = 0.
           END IF
         END DO

#ifdef MPI_PARALL_GRID 
         CALL exchange_p2d(DVDX)
         CALL exchange_p2d(DVDY)
#endif

         IF (.FALSE.) THEN
           OPEN(2305, FILE  = 'erggrad.bin'  , FORM = 'UNFORMATTED') 
           WRITE(2305) 1.
           WRITE(2305) (DVDX(IP), DVDY(IP), SQRT(DVDY(IP)**2+DVDY(IP)**2), IP = 1, MNP)
         ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DIFFERENTIATE_XYDIR_NEIGHADJ(VAR, DVDX, DVDY)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(IN)  :: VAR(MNP)
      REAL(rkind), INTENT(OUT) :: DVDX(MNP), DVDY(MNP)
      INTEGER       :: IP, IADJ, IP_ADJ
      REAL(rkind)   :: xpdiff, ypdiff, dist, eSum
      REAL(rkind)   :: XPcomp, YPcomp, eXP, eYP, eCoeff
      DO IP=1,NP_RES
        eSum=ZERO
        XPcomp=ZERO
        YPcomp=ZERO
        DO IADJ=1,VERT_DEG(IP)
          IP_ADJ=LIST_ADJ_VERT(IADJ,IP)
          xpdiff=XP(IP_ADJ) - XP(IP)
          ypdiff=YP(IP_ADJ) - YP(IP)
          dist=sqrt(xpdiff*xpdiff+ypdiff*ypdiff)
          eSum=eSum + ONE/dist
          eCoeff=(VAR(IP_ADJ) - VAR(IP))/(dist**3)
          XPcomp = XPcomp + xpdiff*eCoeff
          YPcomp = YPcomp + ypdiff*eCoeff
        END DO
        eXP=XPcomp / eSum
        eYP=YPcomp / eSum
        DVDX(IP)=eXP
        DVDY(IP)=eYP
      END DO
#ifdef MPI_PARALL_GRID 
      CALL exchange_p2d(DVDX)
      CALL exchange_p2d(DVDY)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DIFFERENTIATE_XYDIR_AND_SMOOTH(VAR, DVDX, DVDY, nbIter)
      USE DATAPOOL
      IMPLICIT NONE
      integer, intent(in) :: nbIter
      REAL(rkind), INTENT(IN)  :: VAR(MNP)
      REAL(rkind), INTENT(OUT) :: DVDX(MNP), DVDY(MNP)
      INTEGER           :: NI(3)
      INTEGER           :: IE, I1, I2, I3, IP
      REAL(rkind)            :: DEDY(3),DEDX(3)
      REAL(rkind)            :: DVDXIE, DVDYIE
      REAL(rkind)            :: WEI(MNP)
      REAL(rkind)            :: IE_DY(MNE),IE_DX(MNE)
      REAL(rkind)            :: IE_DY_SMO(MNE),IE_DX_SMO(MNE)
      WEI(:)  = 0.0_rkind
      DO IE = 1, MNE 
        NI = INE(:,IE)
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        WEI(NI) = WEI(NI) + 2.*TRIA(IE)
        DEDX(1) = IEN(1,IE)
        DEDX(2) = IEN(3,IE)
        DEDX(3) = IEN(5,IE)
        DEDY(1) = IEN(2,IE)
        DEDY(2) = IEN(4,IE)
        DEDY(3) = IEN(6,IE)
        DVDXIE  = DOT_PRODUCT( VAR(NI),DEDX)
        DVDYIE  = DOT_PRODUCT( VAR(NI),DEDY)
        IE_DX(IE)=DVDXIE
        IE_DY(IE)=DVDYIE
      END DO
      !
      ! Now smoothing
      !
      CALL SMOOTH_ON_TRIANGLE(IE_DX, IE_DX_SMO, nbIter)
      CALL SMOOTH_ON_TRIANGLE(IE_DY, IE_DY_SMO, nbIter)
      !
      ! Mapping to nodes
      !
      DVDX(:) = 0.0_rkind
      DVDY(:) = 0.0_rkind
      DO IE = 1, MNE 
        NI = INE(:,IE)
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        DVDX(NI) = DVDX(NI) + IE_DX_SMO(IE)
        DVDY(NI) = DVDY(NI) + IE_DY_SMO(IE)
      END DO
      DVDX(:) = DVDX(:)/WEI(:)
      DVDY(:) = DVDY(:)/WEI(:)
      DO IP = 1, MNP
        IF (DEP(IP) .LT. DMIN) THEN
          DVDX(IP) = 0.
          DVDY(IP) = 0.
        END IF
      END DO
#ifdef MPI_PARALL_GRID 
      CALL exchange_p2d(DVDX)
      CALL exchange_p2d(DVDY)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVEKCG(DEPLOC, SIGIN, WN, WVC, WVK, WVCG2)
         USE DATAPOOL
         IMPLICIT NONE
         REAL(rkind), INTENT(IN)  :: DEPLOC, SIGIN
         REAL(rkind), INTENT(OUT) :: WVC, WVK, WVCG2, WN
         REAL(rkind) :: SGDLS , AUX1, AUX2
         REAL(rkind) :: WKDEP
! 
         WVC = 0.
         WVK = 0.
         WVCG2 = 0.
         WN = 0.

         IF (SIGIN .LT. VERYSMALL) THEN
            WN = 0.
            WVK=10.
            WVCG2=0.
            RETURN
         END IF

         IF (DEPLOC .GT. DMIN) THEN
            SGDLS = SIGIN*SIGIN*DEPLOC/G9
            AUX1 = 1.0+0.6522*SGDLS+0.4622*(SGDLS**2.0)+0.0864*(SGDLS**4.0)+0.0675*(SGDLS**5.0)
            AUX2 = 1.0/(SGDLS+1.0/AUX1)
            WVC = SQRT(AUX2*G9*DEPLOC)
            WVK = SIGIN/WVC
            WKDEP = WVK*DEPLOC
            IF (WKDEP > 13.0_rkind) THEN
               WN = 0.5_rkind
            ELSE
               WN = 0.5_rkind*(ONE+TWO*WKDEP/SINH(MIN(TWO*KDMAX,TWO*WKDEP)))
            END IF
            WVCG2 = WN*WVC
          ELSE
            WVC  = 0.
            WVK  = 10.
            WVCG2 = 0.
          END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MAKE_WAVE_TABLE()
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IS
         REAL(rkind)    :: SIGIN, WVN, WVC, WVK, WVCG
         REAL(rkind)    :: DEPTH, SIGMAMAX

         NMAX     = IDISPTAB - 1
         DEPTH    = 1.
         SIGMAMAX = SQRT (G9 * DEPFAC)
         DSIGTAB  = SIGMAMAX / MyREAL(NMAX)

         TABK(0)  = 0.
         TABCG(0) = SQRT(G9)

         DO IS = 1, NMAX
           SIGIN = MyREAL(IS)*DSIGTAB
           CALL WAVEKCG(DEPTH, SIGIN, WVN, WVC, WVK, WVCG)
           TABK(IS)  = WVK
           TABCG(IS) = WVCG
         END DO

         IS      = NMAX + 1
         SIGIN   = MyREAL(IS)*DSIGTAB
         CALL WAVEKCG(DEPTH, SIGIN, WVN, WVC, WVK, WVCG)
         TABK(IS)  = WVK
         TABCG(IS) = WVCG

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ALL_FROM_TABLE(SIGIN,DEPTH,WVK,WVCG,WVKDEP,WVN,WVC)

         USE DATAPOOL
         IMPLICIT NONE

         REAL(rkind), INTENT(IN)    :: SIGIN, DEPTH
         REAL(rkind), INTENT(OUT)   :: WVK, WVCG, WVKDEP, WVN, WVC

         REAL(rkind)                :: SQRTDEP, SIGSQDEP, R1, R2, DEPLOC
         INTEGER                    :: ITAB, ITAB2

         DEPLOC   = MAX(DMIN,DEPTH)
         SQRTDEP  = MySQRT(DEPLOC)
         SIGSQDEP = SIGIN * SQRTDEP 
         ITAB     = INT(SIGSQDEP/DSIGTAB)
!
         IF (ITAB.LE.NMAX.AND.ITAB.GE.0) THEN
           ITAB2   = ITAB + 1
           R1      = SIGSQDEP/DSIGTAB - MyREAL(ITAB)
           R2      = ONE - R1
           WVK     = ( R2*TABK(ITAB)  + R1*TABK(ITAB2) ) / DEPLOC 
           WVCG    = ( R2*TABCG(ITAB) + R1*TABCG(ITAB2) ) * SQRTDEP 
           WVKDEP  = WVK * DEPLOC 
           WVN     = ZEROFIVE*(ONE+TWO*WVKDEP/MySINH(MIN(KDMAX,TWO*WVKDEP)))
           WVC     = WVCG/WVN
         ELSE
           WVK     = SIGIN*SIGIN/G9
           WVCG    = ZEROFIVE*G9/SIGIN
           WVKDEP  = WVK * DEPLOC 
           WVN     = ZEROFIVE 
           WVC     = TWO*WVCG
         END IF
!
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CG_FROM_TABLE(SIGIN,DEPTH,WVCG)

         USE DATAPOOL
         IMPLICIT NONE

         REAL(rkind), INTENT(IN)    :: SIGIN, DEPTH
         REAL(rkind), INTENT(OUT)   :: WVCG

         REAL(rkind)                :: SQRTDEP, SIGSQDEP, R1, DEPLOC
         INTEGER             :: ITAB

         DEPLOC   = MAX(DMIN,DEPTH)
         SQRTDEP  = SQRT(DEPLOC)
         SIGSQDEP = SIGIN*SQRTDEP
         ITAB     = INT(SIGSQDEP/DSIGTAB)
         IF (ITAB.LE.NMAX.AND.ITAB.GE.0) THEN
           R1      = SIGSQDEP/DSIGTAB - MyREAL(ITAB)
           WVCG    = ((1.-R1)*TABCG(ITAB)+R1*TABCG(ITAB+1))*SQRTDEP
         ELSE
           WVCG    = .5*G9/SIGIN
         END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE K_FROM_TABLE(SIGIN,DEPTH,WVK)

         USE DATAPOOL
         IMPLICIT NONE

         REAL(rkind), INTENT(IN)    :: SIGIN, DEPTH
         REAL(rkind), INTENT(OUT)   :: WVK

         REAL(rkind)                :: SQRTDEP, SIGSQDEP, R1, DEPLOC
         INTEGER             :: ITAB

         DEPLOC   = MAX(DMIN,DEPTH)
         SQRTDEP  = SQRT(DEPLOC)
         SIGSQDEP = SIGIN * SQRTDEP
         ITAB     = INT(SIGSQDEP/DSIGTAB)
         IF (ITAB.LE.NMAX.AND.ITAB.GE.0) THEN
           R1      = SIGSQDEP/DSIGTAB - MyREAL(ITAB)
           WVK     = ((1.-R1)*TABK(ITAB)+R1*TABK(ITAB + 1))/DEPLOC
         ELSE
           WVK     = SIGIN*SIGIN/G9
         END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_K_C_CG
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER        :: IP, IS
         REAL(rkind)    :: DEPLOC
         REAL(rkind)    :: WVK,WVCG,WVKDEP,WVN,WVC,SPSIGLOC

         DO IP = 1, MNP
           DEPLOC = MAX(DMIN,DEP(IP))
           DO IS = 1, MSC
             SPSIGLOC = SPSIG(IS)
             CALL ALL_FROM_TABLE(SPSIGLOC,DEPLOC,WVK,WVCG,WVKDEP,WVN,WVC)
!             CALL WAVEKCG(DEPLOC,SPSIGLOC,WVN,WVC,WVK,WVCG)
             WK(IS,IP) = WVK
             CG(IS,IP) = WVCG
             WC(IP,IS) = WVC
           END DO
         END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CHECK_STEADY(TIME, CONV1, CONV2, CONV3, CONV4, CONV5)
        USE DATAPOOL
        IMPLICIT NONE
!
!AR: Joseph please check this code ...
!
        REAL(rkind), INTENT(IN)  :: TIME
        REAL(rkind), INTENT(OUT) :: CONV1, CONV2, CONV3, CONV4, CONV5

        INTEGER :: I, IP, IE, IS, ID, NI(3)
        INTEGER :: IPCONV1, IPCONV2, IPCONV3, IPCONV4, IPCONV5, ISCONV(MNP)
        REAL(rkind)  :: SUMAC, ACLOC(MSC,MDC)
        REAL(rkind)  :: ETOT, EAD, DS, HS2, KD
        REAL(rkind)  :: ETOTF3, ETOTF4, TP, KHS2, EFTOT, TM02
        REAL(rkind)  :: FP, CP, KPP, CGP, WNP, UXD, OMEG, OMEG2
        REAL(rkind)  :: CONVK1, CONVK2, CONVK3, CONVK4, CONVK5
        INTEGER :: ITMP

        IPCONV1 = 0
        IPCONV2 = 0
        IPCONV3 = 0
        IPCONV4 = 0
        IPCONV5 = 0

#ifndef MPI_PARALL_GRID
        DO IP = 1, MNP
#else
        DO IP = 1, NP_RES 

          IF(ASSOCIATED(IPGL(IPLG(IP))%NEXT)) THEN !interface node
            IF(IPGL(IPLG(ip))%NEXT%RANK < MYRANK) CYCLE !already in the sum so skip
          ENDIF 
#endif
          ACLOC(:,:) = AC2(:,:,IP)
          SUMAC = SUM(ACLOC)

          ETOT = 0.0
          DO ID = 1, MDC
            DO IS = 2, MSC
              DS = SPSIG(IS) - SPSIG(IS-1)
              EAD = 0.5*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS*DDIR
              ETOT = ETOT + EAD
            END DO
          END DO
          HS2 = 4.0_rkind*SQRT(ETOT)

          ETOTF3 = 0.
          ETOTF4 = 0.
          DO IS = 1, MSC
            DO ID = 1, MDC
              ETOTF3 = ETOTF3 + SPSIG(IS) * ACLOC(IS,ID)**4 * DDIR * DS_BAND(IS)
              ETOTF4 = ETOTF4 +             ACLOC(IS,ID)**4 * DDIR * DS_BAND(IS)
            END DO
          END DO

          IF(ETOTF4 .GT. THR8 .AND. ETOTF3 .GT. THR8) THEN
             FP   = ETOTF3/ETOTF4
 !            CALL WAVEKCG(DEP(IP), FP, WNP, CP, KPP, CGP)
             CALL ALL_FROM_TABLE(FP,DEP(IP),KPP,CGP,KD,WNP,CP)
             TP   = 1.0_rkind/FP/PI2
             KHS2 = HS2 * KPP
          ELSE
             KHS2 = 0.0_rkind
          END IF

          ETOT  = 0.
          EFTOT = 0.
          DO ID=1, MDC
            IF (LSECU .OR. LSTCU) THEN
              UXD  = CURTXY(IP,1)*COSTH(ID) + CURTXY(IP,2)*SINTH(ID)
            ENDIF
            DO IS=1,MSC
              EAD  = SIGPOW(IS,2) * ACLOC(IS,ID) * FRINTF
              IF (LSECU .OR. LSTCU) THEN
                OMEG  = SPSIG(IS) + WK(IS,IP) * UXD
                OMEG2 = OMEG**2
              ELSE
                OMEG2 = SIGPOW(IS,2)
              ENDIF
              ETOT  = ETOT + EAD
              EFTOT = EFTOT + EAD * OMEG2
            ENDDO
          ENDDO

          IF (ETOT/EFTOT .GT. THR) THEN
             TM02 = PI2 * SQRT(ETOT/EFTOT)
          ELSE
             TM02 = 0.
          END IF

          IF (DEP(IP) .LT. DMIN .OR. SUMAC .LT. THR .OR. IOBP(IP) .EQ. 2 .OR. HS2 .LT. THR) THEN
            IPCONV1 = IPCONV1 + 1 ! Summation of the converged grid points ...
            IPCONV2 = IPCONV2 + 1
            IPCONV3 = IPCONV3 + 1
            IPCONV4 = IPCONV4 + 1
            IPCONV5 = IPCONV5 + 1
            ISCONV(IP) = 1
            !write(dbg%fhndl,*) 'boundary -------1------', TIME, IP, IP_IS_STEADY(IP)
            CYCLE
          ELSE 
            CONVK1 = ABS(HSOLD(IP)-HS2)/HS2
            CONVK2 = ABS(HS2-HSOLD(IP))
            CONVK3 = ABS(SUMACOLD(IP)-SUMAC)/SUMAC
            CONVK4 = ABS(KHS2-KHSOLD(IP))/KHSOLD(IP)
            CONVK5 = ABS(TM02-TM02OLD(IP))/TM02OLD(IP)
            IF (CONVK1 .LT. EPSH1) IPCONV1 = IPCONV1 + 1 
            IF (CONVK2 .LT. EPSH2) IPCONV2 = IPCONV2 + 1
            IF (CONVK3 .LT. EPSH3) IPCONV3 = IPCONV3 + 1
            IF (CONVK4 .LT. EPSH4) IPCONV4 = IPCONV4 + 1
            IF (CONVK5 .LT. EPSH5) IPCONV5 = IPCONV5 + 1
            IF (CONVK1 .LT. EPSH1 .AND. CONVK2 .LT. EPSH2 .AND. CONVK3 .LT. EPSH3 .AND. CONVK4 .LT. EPSH4 .AND. CONVK5 .LT. EPSH5) THEN
              ISCONV(IP) = 1 
              !write(dbg%fhndl,*) 'converged -------2------', TIME, IP, IP_IS_STEADY(IP)
            ELSE
              ISCONV(IP) = 0
              !write(dbg%fhndl,*) 'not converged -------3------', TIME, IP, IP_IS_STEADY(IP)
            ENDIF
          END IF
          HSOLD(IP)    = HS2
          SUMACOLD(IP) = SUMAC
          KHSOLD(IP)   = KHS2
          TM02OLD(IP)  = TM02
        END DO  ! IP

        !write(*,*) time, maxval(IP_IS_STEADY), minval(IP_IS_STEADY)

        DO IE = 1, MNE
          NI = INE(:,IE)
          IF (SUM(ISCONV(NI)) .EQ. 3) THEN
            IE_IS_STEADY(IE) = IE_IS_STEADY(IE) + 1
          ELSE
            IE_IS_STEADY(IE) = 0
          ENDIF
          !WRITE(*,*) IE, IE_IS_STEADY(IE)
        ENDDO

        DO IP = 1, MNP
          ITMP = 0
          DO I = 1, CCON(IP)
            IF (IE_IS_STEADY(IE_CELL2(IP,I)) .GE. 1) ITMP = ITMP + 1
          ENDDO
          IF (ITMP .EQ. CCON(IP)) THEN
            IP_IS_STEADY(IP) = IP_IS_STEADY(IP) + 1
            !IF (IP_IS_STEADY(IP) .GT. 2) WRITE(*,*) TIME, IP, IP_IS_STEADY(IP)
          ELSE
            IP_IS_STEADY(IP) = 0
          ENDIF
        ENDDO

#ifdef MPI_PARALL_GRID
        CALL MPI_ALLREDUCE(IPCONV1, itmp, 1, itype, MPI_SUM, COMM, ierr)
        IPCONV1 = itmp
        CALL MPI_ALLREDUCE(IPCONV2, itmp, 1, itype, MPI_SUM, COMM, ierr)
        IPCONV2 = itmp
        CALL MPI_ALLREDUCE(IPCONV3, itmp, 1, itype, MPI_SUM, COMM, ierr)
        IPCONV3 = itmp
        CALL MPI_ALLREDUCE(IPCONV4, itmp, 1, itype, MPI_SUM, COMM, ierr)
        IPCONV4 = itmp
        CALL MPI_ALLREDUCE(IPCONV5, itmp, 1, itype, MPI_SUM, COMM, ierr)
        IPCONV5 = itmp
        CONV1 = MyREAL(IPCONV1)/MyREAL(NP_GLOBAL)*100.0
        CONV2 = MyREAL(IPCONV2)/MyREAL(NP_GLOBAL)*100.0
        CONV3 = MyREAL(IPCONV3)/MyREAL(NP_GLOBAL)*100.0
        CONV4 = MyREAL(IPCONV4)/MyREAL(NP_GLOBAL)*100.0
        CONV5 = MyREAL(IPCONV5)/MyREAL(NP_GLOBAL)*100.0
#else
        CONV1 = MyREAL(IPCONV1)/MyREAL(MNP)*100.0
        CONV2 = MyREAL(IPCONV2)/MyREAL(MNP)*100.0
        CONV3 = MyREAL(IPCONV3)/MyREAL(MNP)*100.0
        CONV4 = MyREAL(IPCONV4)/MyREAL(MNP)*100.0
        CONV5 = MyREAL(IPCONV5)/MyREAL(MNP)*100.0
#endif

#ifdef MPI_PARALL_GRID
         IF (myrank == 0) THEN
           WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 1 REACHED IN', CONV1, '% GRIDPOINTS'
           WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 2 REACHED IN', CONV2, '% GRIDPOINTS'
           WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 3 REACHED IN', CONV3, '% GRIDPOINTS'
           WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 4 REACHED IN', CONV4, '% GRIDPOINTS'
           WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 5 REACHED IN', CONV5, '% GRIDPOINTS'
         END IF
#else
         WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 1 REACHED IN', CONV1, '% GRIDPOINTS'
         WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 2 REACHED IN', CONV2, '% GRIDPOINTS'
         WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 3 REACHED IN', CONV3, '% GRIDPOINTS'
         WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 4 REACHED IN', CONV4, '% GRIDPOINTS'
         WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 5 REACHED IN', CONV5, '% GRIDPOINTS'
#endif
         FLUSH(STAT%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
     SUBROUTINE TWOD2ONED(ACLOC,AC1D)
         USE DATAPOOL
         IMPLICIT NONE

         REAL(rkind), INTENT(IN)   :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(OUT)  :: AC1D(MSC*MDC)

         INTEGER            :: IS, ID


         DO IS = 1, MSC
           DO ID = 1, MDC
             AC1D(ID + (IS-1) * MDC) = ACLOC(IS,ID)
           END DO
         END DO

     END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
     SUBROUTINE ONED2TWOD(AC1D,ACLOC)
         USE DATAPOOL
         IMPLICIT NONE

         REAL(rkind), INTENT(OUT)   :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(IN)  :: AC1D(MSC*MDC)

         INTEGER            :: IS, ID


         DO IS = 1, MSC
           DO ID = 1, MDC
             ACLOC(IS,ID) = AC1D(ID + (IS-1) * MDC)
           END DO
        END DO

     END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
!      REAL(rkind) FUNCTION GAMMA_FUNC(XX)
      FUNCTION GAMMA_FUNC(XX)
         USE DATAPOOL, ONLY: rkind
         IMPLICIT NONE
         REAL(rkind) :: GAMMA_FUNC
!
!     Purpose:
!        Compute the transcendental function Gamma
!
!     Subroutines used
!        GAMMLN  (Numerical Recipes)
!
         real(rkind) GAMMLN
         REAL(rkind) XX, YY, ABIG
         SAVE ABIG
         DATA ABIG /30./

         YY = GAMMLN(XX)
         IF (YY > ABIG) YY = ABIG
         IF (YY < -ABIG) YY = -ABIG
         GAMMA_FUNC = EXP(YY)

      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      FUNCTION GAMMLN(XX)
         USE DATAPOOL, ONLY: rkind, ONEHALF, ONE
         IMPLICIT NONE
         real(rkind) :: GAMMLN
!
!     Method:
!        function is copied from: Press et al., "Numerical Recipes"
!
         real(rkind) XX
         INTEGER J
         real(rkind)  COF(6),STP,FPF,X,TMP,SER
         DATA COF,STP/76.18009173_rkind,-86.50532033_rkind,       &
     &               24.01409822_rkind,-1.231739516_rkind,        &
     &               .120858003E-2_rkind,-.536382E-5_rkind,       &
     &               2.50662827465_rkind/
         DATA FPF/5.5_rkind/
         X = XX-ONE
         TMP = X+FPF
         TMP = (X+ONEHALF)*LOG(TMP)-TMP
         SER = ONE
         DO J = 1, 6
            X = X+ONE
           SER = SER+COF(J)/X
         END DO
         GAMMLN = TMP+LOG(STP*SER)

      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      FUNCTION VEC2RAD(U,V)
         USE DATAPOOL, ONLY: rkind, pi
         IMPLICIT NONE
         REAL(RKIND) :: VEC2RAD

         REAL(rkind)            :: U,V

         VEC2RAD = MyATAN2(V,U) * 180/PI
         IF (VEC2RAD < 0.0) VEC2RAD = VEC2RAD + 360.0_rkind
         VEC2RAD = VEC2RAD * PI/180.

      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      FUNCTION VEC2DEG(U,V)
         USE DATAPOOL, ONLY: rkind, pi 
         IMPLICIT NONE
         REAL(RKIND) :: VEC2DEG
         REAL(rkind), intent(in) :: U,V

         VEC2DEG = MyATAN2(V,U) * 180./PI
         IF (VEC2DEG < 0.0) VEC2DEG = VEC2DEG + 360.0_rkind

      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      FUNCTION DVEC2RAD(U,V)
         USE DATAPOOL, ONLY: rkind, pi
         IMPLICIT NONE
         REAL(rkind)            :: U,V
         REAL(rkind)  :: DVEC2RAD

         DVEC2RAD = MyATAN2(V,U) * 180.0_rkind/PI
         IF (DVEC2RAD < 0.0_rkind) DVEC2RAD = DVEC2RAD + 360.0_rkind
         DVEC2RAD = DVEC2RAD * PI/180.0_rkind

      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      FUNCTION DVEC2DEG(U,V)
         USE DATAPOOL, ONLY : PI, rkind
         IMPLICIT NONE
         REAL(rkind) :: DVEC2DEG

         REAL(rkind)            :: U,V

         DVEC2DEG = MyATAN2(V,U) * 180.0_rkind/PI
         IF (DVEC2DEG < 0.0_rkind) DVEC2DEG = DVEC2DEG + 360.0_rkind

      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEG2NAUT (DEGREE, DEG, LTRANS)
         USE DATAPOOL, ONLY: rkind
         IMPLICIT NONE
!
!           Nautical convention           Cartesian convention
!     (Where the wind/waves come from)  (Where the wind/waves go to)
!                    0                             90
!                    |                              |
!                    |                              |
!                    |                              |
!                    |                              |
!        270 --------+-------- 90       180 --------+-------- 0
!                    |                              |
!                    |                              |
!                    |                              |
!                    |                              |
!                   180                            270
!
      LOGICAL  :: LTRANS
      REAL(rkind), INTENT(IN)     :: DEGREE
      REAL(rkind), INTENT(OUT)    :: DEG

      REAL(rkind)                 :: DNORTH

      IF ( LTRANS ) THEN
          DNORTH = 90.0
          DEG    = 180. + DNORTH - DEGREE
      ELSE
          DEG    = DEGREE
      END IF
!
      IF (DEG .GE. 360.0_rkind) THEN
        DEG = MOD (DEG, 360.0_rkind)
      ELSE IF (DEG .LT. 0.) THEN
        DEG = MOD (DEG, 360.0_rkind) + 360.0_rkind
      ELSE
!       DEG between 0 and 360; do nothing
      endif
!
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       logical function my_isnan(a)
         use datapool, only : rkind
         real(rkind) :: a
         if (a.ne.a) then
           my_isnan = .true.
         else
           my_isnan = .false.
         end if
         return
         end
!**********************************************************************
!*                                                                    *
!**********************************************************************
         logical function my_isinf(a)
         use datapool, only : rkind
         real(rkind) :: a
         if (int(a*0) .ne. 0) then
           my_isinf = .true.
         else
           my_isinf = .false.
         end if
         return
         end
!**********************************************************************
!*                                                                    *
!**********************************************************************
         logical function iseq0(a)
         use datapool, only : rkind
         real(rkind) :: a

         if (abs(a) .lt. tiny(1.0)) then
           iseq0 = .true.
         else
           iseq0 = .false.
         end if
         return
         end
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CHKVOL(VAL,PO,NE,POSNEG)
         USE DATAPOOL
         IMPLICIT NONE

         REAL(rkind), INTENT(IN)    :: VAL(MNP)
         REAL(rkind), INTENT(OUT)   :: PO,NE,POSNEG

         REAL(rkind)                :: TMP
         INTEGER               :: IE,I1,I2,I3

         PO = ZERO
         NE = ZERO
         DO IE = 1, MNE
            I1 = INE(1,IE)
            I2 = INE(2,IE)
            I3 = INE(3,IE)
            TMP = ONETHIRD * (VAL(I1)+VAL(I2)+VAL(I3)) * TRIA(IE)
            IF (TMP .GT. ZERO) THEN
              PO = PO + TMP
            ELSE
              NE = NE + TMP
            END IF
            POSNEG = PO + NE
         END DO

       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CHECKCONS(VAL,SUMAC)

         USE DATAPOOL
         IMPLICIT NONE

         REAL(rkind), INTENT(IN)    :: VAL(MNP)
         REAL(rkind), INTENT(OUT)   :: SUMAC

         REAL(rkind)                :: TMP
         INTEGER               :: IE,I1,I2,I3

         SUMAC = 0.

         DO IE = 1, MNE

            I1 = INE(1,IE)
            I2 = INE(2,IE)
            I3 = INE(3,IE)

            TMP = 1./3. * (VAL(I1)+VAL(I2)+VAL(I3)) * TRIA(IE)

            SUMAC = SUMAC + TMP

         END DO

       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ADVTEST(INIT)
        USE DATAPOOL
        IMPLICIT NONE

        INTEGER :: IP
        REAL(rkind)  :: R
        REAL(rkind), INTENT(OUT)   :: INIT(MNP)

        INIT = 0.0_rkind
        DO IP = 1, MNP
          R = SQRT( (XP(IP)+ONEHALF)**2 + YP(IP)**2 )
          IF ( R <= 0.25_rkind ) THEN
            IF (LZYLINDER) THEN
              INIT(IP) = ONE
            ELSE
              INIT(IP) = COS(PI2*R)
            END IF
          ELSE
            INIT(IP) = ZERO
          END IF
        END DO

        WRITE(4001)  SNGL(RTIME)
        WRITE(4001) (1., 1., SNGL(INIT(IP)), IP = 1, MNP)

       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ERG2WWM(STEPS)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER :: T, IP, STEPS
         REAL(rkind)    :: HEADSP
         REAL(rkind), PARAMETER :: FRFAK = 0.9
         REAL(rkind)    :: VEC2RAD, ANG, VEL, WAVEL, DEPTH
         REAL(rkind), ALLOCATABLE  :: HP(:), QU(:), QV(:), UP(:), VP(:)

         HEADSP = 0.0

#ifdef MPI_PARALL_GRID
         CALL WWM_ABORT('ERG2WWM CANNOT BE CALLED FROM SCHISM')
#endif

         ALLOCATE (QU(MNP), QV(MNP), UP(MNP), VP(MNP), HP(MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_aux, allocate error 1')

         QU = 0.
         QV = 0.
         UP = 0.
         VP = 0.
         HP = 0.

         OPEN(1230, FILE = 'quqvh.bin', FORM = 'UNFORMATTED')
         OPEN(1240, FILE = 'current.dat')
         OPEN(1250, FILE = 'wlevel.dat')

         DO T = 1, STEPS

           READ(1230) HEADSP
           READ(1230) (QU(IP), QV(IP), HP(IP) , IP = 1, MNP)
           WRITE (1240, *) HEADSP
           WRITE (1250, *) HEADSP

           DO IP = 1, MNP
             DEPTH = DEP(IP)+HP(IP)
             IF ( DEPTH .GT. DMIN ) THEN
               UP(IP) = QU(IP)/DEPTH
               VP(IP) = QV(IP)/DEPTH
               VEL = SQRT(UP(IP)**2 + VP(IP)**2)
               ANG = VEC2RAD(UP(IP),VP(IP))
               WAVEL = SQRT(G9*DEPTH)
               IF (VEL .GT. FRFAK * WAVEL) THEN
                 UP(IP) = COS(ANG)*FRFAK*WAVEL
                 VP(IP) = SIN(ANG)*FRFAK*WAVEL
               END IF
             ELSE
               UP(IP) = 0.0
               VP(IP) = 0.0
             END IF
           END DO

           WRITE (1240, 1200) UP(:)
           WRITE (1240, 1200) VP(:)
           WRITE (1250, 1200) HP(:)

         END DO

         DEALLOCATE (QU,QV,HP,UP,VP)

1200     FORMAT(8F9.4)

         CLOSE (1240)
         CLOSE (1250)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTELEMENT_COEF(X,Y,XP,YP,WI)
      USE DATAPOOL, ONLY : RKIND
      IMPLICIT NONE
      REAL(rkind),    INTENT(IN)  :: X(3), Y(3)
      REAL(rkind),    INTENT(IN)  :: XP, YP
      REAL(rkind),  INTENT(OUT) :: WI(3)

      REAL(rkind) :: y1,y2,y3,x1,x2,x3
      REAL(rkind) :: n1, d1, n2, d2, n3, d3
      x1 = X(1)
      x2 = X(2)
      x3 = X(3)
      y1 = Y(1)
      y2 = Y(2)
      y3 = Y(3)
      n1=(XP-x2)*(y3-y2) - (YP-y2)*(x3-x2)
      d1=(x1-x2)*(y3-y2) - (y1-y2)*(x3-x2)
      n2=(XP-x1)*(y3-y1) - (YP-y1)*(x3-x1)
      d2=(x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)
      n3=(XP-x1)*(y2-y1) - (YP-y1)*(x2-x1)
      d3=(x3-x1)*(y2-y1) - (y3-y1)*(x2-x1)
      Wi(1)=n1/d1
      Wi(2)=n2/d2
      Wi(3)=n3/d3
      END SUBROUTINE INTELEMENT_COEF
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTELEMENT(X,Y,Z,XP,YP,Wi,Zi,LSAME)
      USE DATAPOOL, ONLY : RKIND, THR8
      IMPLICIT NONE

      LOGICAL, INTENT(IN)  :: LSAME
      REAL(rkind),    INTENT(IN)  :: X(3), Y(3), Z(3)
      REAL(rkind),    INTENT(IN)  :: XP, YP
      REAL(rkind),    INTENT(OUT) :: Zi
      REAL(rkind),  INTENT(OUT) :: WI(3)

      REAL(rkind)               :: y1,y2,y3,x1,x2,x3,z1,z2,z3
      REAL(rkind), SAVE         :: A,B,C,D

      IF (.NOT. LSAME) THEN
        x1 = X(1)
        x2 = X(2)
        x3 = X(3)
        y1 = Y(1)
        y2 = Y(2)
        y3 = Y(3)
        z1 = Z(1)
        z2 = Z(2)
        z3 = Z(3)
        A = y1*(z2 - z3)  +  y2*(z3 - z1) +  y3*(z1 - z2)
        B = z1*(x2 - x3)  +  z2*(x3 - x1) +  z3*(x1 - x2)
        C = x1*(y2 - y3)  +  x2*(y3 - y1) +  x3*(y1 - y2)
        D = -A*x1 - B*y1 - C*z1
        IF (ABS(C) .GT. THR8 ) THEN
          WI(1) = -A/C
          WI(2) = -B/C
          WI(3) = -D/C
        ELSE
          WI    = 0.0_rkind
        endif
      END IF
      Zi = WI(1) * XP + WI(2) * YP + WI(3)
 
      END SUBROUTINE INTELEMENT
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTELEMENT_AC_LOC(I,ACLOC,CURTXYLOC,DEPLOC,WATLEVLOC,WKLOC)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: I
      REAL(rkind), INTENT(OUT) :: ACLOC(MSC,MDC)
      REAL(rkind), INTENT(OUT) :: CURTXYLOC(2), DEPLOC, WATLEVLOC, WKLOC(MSC)
      REAL(rkind), SAVE         :: WI(3)
      INTEGER :: IS, NI(3), IE
      REAL(rkind) :: WVN, WVC, WVK, WVCG, WVKDEP
      integer IP, J

      IE = STATION(I)%ELEMENT
      WI = STATION(I)%WI
      NI = INE(:,IE)

      DEPLOC    = ZERO 
      WATLEVLOC = ZERO 
      CURTXYLOC = ZERO 
      ACLOC     = ZERO 

      DO J=1,3
        IP           = INE(J,IE)
        DEPLOC       = DEPLOC       + WI(J) * DEP(IP)
        CURTXYLOC(:) = CURTXYLOC(:) + WI(J) * CURTXY(IP,:)
        WATLEVLOC    = WATLEVLOC    + WI(J) * WATLEV(IP)
        ACLOC(:,:)   = ACLOC(:,:)   + WI(J) * AC2(:,:,IP)
        !write(*,'(3I10,5F15.8)') j, ie, ip, wi(j), DEP(IP), WATLEV(IP), CURTXY(IP,:) 
      END DO

      DO IS = 1, MSC
        !CALL WAVEKCG(DEPLOC, SPSIG(IS), WVN, WVC, WVK, WVCG)
        CALL ALL_FROM_TABLE(SPSIG(IS),DEPLOC,WVK,WVCG,WVKDEP,WVN,WVC)
        WKLOC(IS) = WVK
      END DO

      END SUBROUTINE INTELEMENT_AC_LOC
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTELEMENT_WW3GLOBAL_LOC(IE,XPC,YPC,WW3LOCAL)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: IE
      REAL(rkind), INTENT(IN)  :: XPC, YPC
      REAL(rkind), INTENT(OUT) :: WW3LOCAL(8)
      REAL(rkind)          :: XOUTELE(3), YOUTELE(3)
      REAL(rkind), SAVE    :: WI(3)
      INTEGER     :: I, NI(3)
      LOGICAL :: LSAME

      NI       = INE(:,IE)
      XOUTELE  = XP(NI)
      YOUTELE  = YP(NI)

      LSAME = .FALSE.

      DO I = 1, 8 
         CALL INTELEMENT(XOUTELE,YOUTELE,WW3GLOBAL(I,NI),XPC,YPC,WI,WW3LOCAL(I),LSAME)
         LSAME = .TRUE.
      END DO

      END SUBROUTINE INTELEMENT_WW3GLOBAL_LOC
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CSEVAL ( NFU, FILEN, LSE, DIMS, SEVAL, MULTIPLE_IN)
      USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)          :: NFU
         LOGICAL, INTENT(IN)          :: MULTIPLE_IN
         CHARACTER(LEN=*), INTENT(IN) :: FILEN

         LOGICAL, INTENT(IN)          :: LSE

         INTEGER, INTENT(IN)          :: DIMS

         REAL(rkind), INTENT(INOUT)   :: SEVAL(MNP, DIMS)
         REAL(rkind)                  :: SEVAL2(NP_TOTAL, DIMS)
#ifdef MPI_PARALL_GRID
         INTEGER                      :: IP
         REAL(rkind)                  :: Vtotal(np_total), Vlocal(MNP)
#endif
         INTEGER                      :: IC, IFSTAT, IDIM
         CHARACTER(LEN=128)           :: HEADLN
#ifdef MPI_PARALL_GRID
         IF (MULTIPLE_IN .or. (myrank .eq. 0)) THEN
#endif
           IF (LSE) THEN
             READ(NFU,*) HEADLN
             WRITE(STAT%FHNDL,'("+TRACE...",2A)') 'Reading the header of the serial file ... HEADER ', TRIM(HEADLN)
             DO IC = 1, DIMS
               READ( NFU, *, IOSTAT = IFSTAT ) SEVAL2(:, IC)
               IF ( IFSTAT /= 0 ) CALL WWM_ABORT('unexpected error reading the serial file in CSEVAL 1')
             END DO
           ELSE
             OPEN( NFU, FILE = TRIM(FILEN), STATUS = 'OLD')
             WRITE(STAT%FHNDL,'("+TRACE...",2A)') 'Reading the file of the request ... ', TRIM(FILEN)
             READ(NFU,*) HEADLN
             DO IC = 1, DIMS
               READ( NFU, *, IOSTAT = IFSTAT ) SEVAL2(:, IC)
               IF ( IFSTAT /= 0 ) CALL WWM_ABORT(' unexpected error reading the serial file in CSEVAL 2')
             END DO
             CLOSE( NFU )
           END IF
#ifdef MPI_PARALL_GRID
         END IF
#endif
#ifdef MPI_PARALL_GRID
         IF (MULTIPLE_IN) THEN
           DO IP=1,MNP
             SEVAL(IP,:) = SEVAL2(IPLG(IP),:)
           END DO
         ELSE
           DO IDIM=1,DIMS
             IF (myrank .eq. 0) THEN
               Vtotal=SEVAL2(:,IDIM)
             END IF
             CALL SCATTER_ONED_ARRAY(Vtotal, Vlocal)
             SEVAL(:,IDIM)=Vlocal
           END DO
         END IF
#else
         SEVAL(:,:) = SEVAL2(:,:)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
        SUBROUTINE ERROR(X,ERR)
        USE DATAPOOL, ONLY : rkind,pi
        IMPLICIT NONE
        REAL(rkind) :: EPS,X,X2,ERR,ER,R,C0
        INTEGER     :: K
!       =========================================
!       Purpose: Compute error function erf(x)
!       Input:   x   --- Argument of erf(x)
!       Output:  ERR --- erf(x)
!       =========================================
        EPS=1.E-15_rkind
        X2=X*X
        IF (MyABS(X).LT.3.5D0) THEN
           ER=1.0D0
           R=1.0D0
           DO 10 K=1,50
              R=R*X2/(K+0.5D0)
              ER=ER+R
              IF (MyABS(R).LE.MyABS(ER)*EPS) GO TO 15
10         CONTINUE
15         C0=2.0D0/MySQRT(PI)*X*MyEXP(-X2)
           ERR=C0*ER
        ELSE
           ER=1.0D0
           R=1.0D0
           DO K=1,12
              R=-R*(K-0.5D0)/X2
              ER=ER+R
           END DO
           C0=MyEXP(-X2)/(MyABS(X)*MySQRT(PI))
           ERR=1.0D0-C0*ER
           IF (X.LT.0.0) ERR=-ERR
        endif
        END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
        SUBROUTINE WRINPGRD_XFN()
#ifdef MPI_PARALL_GRID
         USE DATAPOOL, ONLY : myrank, comm, ierr, MNP, XP, YP, DEP, MNE, INE
#else
         USE DATAPOOL, ONLY : MNP, XP, YP, DEP, MNE, INE
#endif
         IMPLICIT NONE
         INTEGER :: I

#ifdef MPI_PARALL_GRID
         if (myrank == 0) then
#endif
         OPEN(2222, FILE='systest.dat', STATUS='UNKNOWN')

         CALL XFNHEADER_1(2222,0,MNP)

         DO I = 1, MNP
           WRITE(2222,'(I10,2F20.8,F15.4)') I-1, XP(I), YP(I), DEP(I)
         END DO

         CALL XFNHEADER_2(2222,MNE)

         DO I = 1, MNE
           WRITE(2222,'(5I10)') INE(1,I)-1, INE(2,I)-1, INE(3,I)-1, 0, I-1
         END DO
#ifdef MPI_PARALL_GRID
        endif
#endif

        END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
        SUBROUTINE WRINPGRD_SHP()
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: I

#ifdef MPI_PARALL_GRID
         if (myrank == 0) then
#endif
         OPEN(2222, FILE='systest.dat', STATUS='UNKNOWN')

         CALL XFNHEADER_1(2222,0,MNP)

         DO I = 1, MNP
           WRITE(2222,'(I10,2F20.8,F15.4)') I-1, XP(I), YP(I), DEP(I)
         END DO

         CALL XFNHEADER_2(2222,MNE)

         DO I = 1, MNE
           WRITE(2222,'(5I10)') INE(1,I)-1, INE(2,I)-1, INE(3,I)-1, 0, I-1
         END DO
#ifdef MPI_PARALL_GRID
        endif
#endif

        END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE XFNHEADER_1(IFILE,NKR,NKG)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IFILE, NKR, NKG

      WRITE(IFILE,'(A)') 'C system.dat, made by tri2sys'
      WRITE(IFILE,'(A)') 'C Number of Boundary Nodes:'
      WRITE(IFILE,'(I10)') NKR
      WRITE(IFILE,'(A)') 'C Number of Domain Nodes:'
      WRITE(IFILE,'(I10)') NKG
      WRITE(IFILE,'(A)') 'C Koordinaten und Skalarwerte der Knoten'
      WRITE(IFILE,'(A)') 'C --------------------------------------'
      WRITE(IFILE,'(A)') 'C Zuerst die Randknoten  (Anzahl s.o.),'
      WRITE(IFILE,'(A)') 'C dann die Gebietsknoten (Anzahl s.o.).'
      WRITE(IFILE,'(A)') 'C ------------+-------------+-------------+---------------'
      WRITE(IFILE,'(A)') 'C     Nr.     |  x-Koord.   |   y-Koord.  | Skalarwert'
      WRITE(IFILE,'(A)') 'C ------------+-------------+-------------+---------------'

      END SUBROUTINE XFNHEADER_1
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE XFNHEADER_2(IFILE,NELEM)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IFILE, NELEM

      WRITE(IFILE,'(A)') "C ------------------------------------------------------------"
      WRITE(IFILE,'(A)') "C Anzahl der Elemente:"
      WRITE(IFILE,'(I11)') NELEM
      WRITE(IFILE,'(A)') "C Elementverzeichnis"
      WRITE(IFILE,'(A)') "C ------------------------------------------------------------"
      WRITE(IFILE,'(A)') "C    Knoten i  Knoten j  Knoten k   Kennung     Nr."

      END SUBROUTINE XFNHEADER_2
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE RHEADER_NODE(IFILE,NKR,NKG)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IFILE
      INTEGER, INTENT(OUT):: NKR, NKG

      READ(IFILE,*) 
      READ(IFILE,*)
      READ(IFILE,*) NKR
      READ(IFILE,*)
      READ(IFILE,*) NKG
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)

      END SUBROUTINE RHEADER_NODE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE RHEADER_ELEMENT(IFILE,NELEM)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IFILE
      INTEGER, INTENT(OUT):: NELEM

      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*) NELEM
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)

      END SUBROUTINE RHEADER_ELEMENT
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE FIND_ELE ( M,N,Lelem,Xkno,Xp,Yp,Ele )
      USE DATAPOOL, ONLY : THR8, RKIND
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER, INTENT(IN)      :: M, N
      REAL(rkind), INTENT(IN)  :: Xkno(2,N), Xp , Yp
      INTEGER, INTENT(IN)      :: Lelem(3,M)
      INTEGER, INTENT(INOUT)   :: Ele
!
! Local variables
!
      REAL(rkind), SAVE   :: xi, xj, xk, yi, yj, yk, dx, dy, f, Xp8, Yp8
      REAL(rkind), SAVE   :: xmax , xmin , ymax , ymin
      INTEGER, SAVE  :: i , i0 , idx , ielem , if0 , if1 , ijk , k , ki , kj , kk , l

      DATA idx/0/ , i0/0/ , ielem/1/ , if0/0/ , if1/0/
!
!     Laengste Kannte (DX,DY) bestimmen
!     Dieser Programmabschnitt wird nur beim ersten Aufruf von PLO160
!     durchlaufen!
!
      Xp8 = Xp
      Yp8 = Yp

      IF ( idx/=1 ) THEN
         idx = 1
         DO i = 1 , M
            ki = Lelem(1,i)! + 1
            kj = Lelem(2,i)! + 1
            kk = Lelem(3,i)! + 1
            xi = Xkno(1,ki)
            yi = Xkno(2,ki)
            xj = Xkno(1,kj)
            yj = Xkno(2,kj)
            xk = Xkno(1,kk)
            yk = Xkno(2,kk)
            IF ( i==1 ) THEN
               dx = MAX(ABS(xi-xj),ABS(xi-xk),ABS(xj-xk))
               dy = MAX(ABS(yi-yj),ABS(yi-yk),ABS(yj-yk))
            ELSE
               dx = MAX(dx,ABS(xi-xj),ABS(xi-xk),ABS(xj-xk))
               dy = MAX(dy,ABS(yi-yj),ABS(yi-yk),ABS(yj-yk))
            endif
         ENDDO
      endif
!     ------------------------------------------------------------------
!     TEST, OB DER PUNKT IM ZULETZT ANGESPROCHENEN ELEMENT LIEGT
!     ------------------------------------------------------------------
      IF ( i0==1 .AND. Ele/=0 ) THEN
         IF ( Yp8-ymin > THR8 ) THEN
            IF ( Yp8-ymax < -THR8 ) THEN
               IF ( Xp8-xmin > THR8 ) THEN
                  IF ( Xp8-xmax < -THR8 ) THEN
                     f = xi*(yj-Yp8) + xj*(Yp8-yi) + Xp8*(yi-yj)
                     IF ( f > THR8 ) THEN
                        f = xj*(yk-Yp8) + xk*(Yp8-yj) + Xp8*(yj-yk)
                        IF ( f > THR8  ) THEN
                           f = xk*(yi-Yp8) + xi*(Yp8-yk) + Xp8*(yk-yi)
                           IF ( f > THR8 ) THEN
                              Ele = ielem ! Element gefunden -->RETURN
                              RETURN
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
      endif
!     ------------------------------------------------------------------
!     Element suchen
!     ------------------------------------------------------------------
      i0 = 1
      i = ielem
      IF ( i<1 ) i = 1
      k = i
      l = i
      ijk = 0

100   DO
         ijk = ijk + 1
!.....   ABFRAGE AUF X-Richtung
         ki = Lelem(1,i)! + 1
         xi = Xkno(1,ki)
         IF ( MyABS(xi-Xp8)<=dx ) THEN
            kj = Lelem(2,i)! + 1
            kk = Lelem(3,i)! + 1
            xj = Xkno(1,kj)
            xk = Xkno(1,kk)
!.....    Punkt ausserhalb Element:
            xmin = MIN(xi,xj,xk)
            IF ( Xp8>=xmin ) THEN
               xmax = MAX(xi,xj,xk)
               IF ( Xp8<=xmax ) THEN
!.....        ABFRAGE AUF Y-Richtung
                  yi = Xkno(2,ki)
                  IF ( MyABS(yi-Yp8)<=dy ) THEN
                     yj = Xkno(2,kj)
                     yk = Xkno(2,kk)
!.....          Punkt ausserhalb Element:
                     ymin = MIN(yi,yj,yk)
                     IF ( Yp8>=ymin ) THEN
                        ymax = MAX(yi,yj,yk)
                        IF ( Yp8<=ymax ) THEN
!.....              Bis jetzt liegt Punkt innerhalb des das Element
!                   umschlieszenden Rechtecks XMIN/XMAX, YMIN/YMAX
!                   Pruefen, ob Punkt wirklich innerhalb DREIECK-Element
!                   liegt: BERECHNUNG DER TEILFLAECHEN (ohne 0.5)
                           f = xi*(yj-Yp8) + xj*(Yp8-yi) + Xp8*(yi-yj)
                           IF ( f>=0.0_rkind) THEN
                              f = xj*(yk-Yp8) + xk*(Yp8-yj) + Xp8*(yj-yk)
                              IF ( f>=0.0_rkind ) THEN
                                 f = xk*(yi-Yp8) + xi*(Yp8-yk) + Xp8*(yk-yi)
                                 IF ( f>=0.0_rkind ) THEN
                                    Ele = i
                                    ielem = Ele
                                    RETURN
                                 endif
                              endif
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
!     SCHLEIFE UEBER ALLE ELEMENTE wird hier folgendermassen hochgezaehlt:
!     beginnend bei IEALT, im Wechsel nach vorn und rueckwaerts suchend
         IF ( k<M .AND. if1==0 ) THEN
            if0 = 0
            IF ( l>1 ) if1 = 1
            k = k + 1
            i = k
            IF ( ijk<=M ) CYCLE
         endif
         CONTINUE
         EXIT
      ENDDO

      IF ( l>1 .AND. if0==0 ) THEN
         if1 = 0
         IF ( k<M ) if0 = 1
         l = l - 1
         i = l
      endif

      IF ( ijk<=M ) GOTO 100

      Ele = 0 
 
      END SUBROUTINE FIND_ELE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE FIND_ELE_WIND ( M,N,Lelem,Xkno,Xp,Yp,Ele )
      USE DATAPOOL, ONLY : THR8, RKIND
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER, INTENT(IN)      :: M, N
      REAL(rkind), INTENT(IN)  :: Xkno(2,N), Xp , Yp
      INTEGER, INTENT(IN)      :: Lelem(3,M)
      INTEGER, INTENT(INOUT)   :: Ele
!
! Local variables
!
      REAL(rkind), SAVE   :: xi, xj, xk, yi, yj, yk, dx, dy, f, Xp8, Yp8
      REAL(rkind), SAVE   :: xmax , xmin , ymax , ymin
      INTEGER, SAVE  :: i , i0 , idx , ielem , if0 , if1 , ijk , k , ki , kj , kk , l

      DATA idx/0/ , i0/0/ , ielem/1/ , if0/0/ , if1/0/
!
!     Laengste Kannte (DX,DY) bestimmen
!     Dieser Programmabschnitt wird nur beim ersten Aufruf von PLO160
!     durchlaufen!
!
      Xp8 = Xp
      Yp8 = Yp

      IF ( idx/=1 ) THEN
         idx = 1
         DO i = 1 , M
            ki = Lelem(1,i)! + 1
            kj = Lelem(2,i)! + 1
            kk = Lelem(3,i)! + 1
            xi = Xkno(1,ki)
            yi = Xkno(2,ki)
            xj = Xkno(1,kj)
            yj = Xkno(2,kj)
            xk = Xkno(1,kk)
            yk = Xkno(2,kk)
            IF ( i==1 ) THEN
               dx = MAX(ABS(xi-xj),ABS(xi-xk),ABS(xj-xk))
               dy = MAX(ABS(yi-yj),ABS(yi-yk),ABS(yj-yk))
            ELSE
               dx = MAX(dx,ABS(xi-xj),ABS(xi-xk),ABS(xj-xk))
               dy = MAX(dy,ABS(yi-yj),ABS(yi-yk),ABS(yj-yk))
            endif
         ENDDO
      endif
!     ------------------------------------------------------------------
!     TEST, OB DER PUNKT IM ZULETZT ANGESPROCHENEN ELEMENT LIEGT
!     ------------------------------------------------------------------
      IF ( i0==1 .AND. Ele/=0 ) THEN
         IF ( Yp8-ymin > THR8 ) THEN
            IF ( Yp8-ymax < -THR8 ) THEN
               IF ( Xp8-xmin > THR8 ) THEN
                  IF ( Xp8-xmax < -THR8 ) THEN
                     f = xi*(yj-Yp8) + xj*(Yp8-yi) + Xp8*(yi-yj)
                     IF ( f > THR8 ) THEN
                        f = xj*(yk-Yp8) + xk*(Yp8-yj) + Xp8*(yj-yk)
                        IF ( f > THR8  ) THEN
                           f = xk*(yi-Yp8) + xi*(Yp8-yk) + Xp8*(yk-yi)
                           IF ( f > THR8 ) THEN
                              Ele = ielem ! Element gefunden -->RETURN
                              RETURN
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
      endif
!     ------------------------------------------------------------------
!     Element suchen
!     ------------------------------------------------------------------
      i0 = 1
      i = ielem
      IF ( i<1 ) i = 1
      k = i
      l = i
      ijk = 0

100   DO
         ijk = ijk + 1
!.....   ABFRAGE AUF X-Richtung
         ki = Lelem(1,i)! + 1
         xi = Xkno(1,ki)
         IF ( MyABS(xi-Xp8)<=dx ) THEN
            kj = Lelem(2,i)! + 1
            kk = Lelem(3,i)! + 1
            xj = Xkno(1,kj)
            xk = Xkno(1,kk)
!.....    Punkt ausserhalb Element:
            xmin = MIN(xi,xj,xk)
            IF ( Xp8>=xmin ) THEN
               xmax = MAX(xi,xj,xk)
               IF ( Xp8<=xmax ) THEN
!.....        ABFRAGE AUF Y-Richtung
                  yi = Xkno(2,ki)
                  IF ( MyABS(yi-Yp8)<=dy ) THEN
                     yj = Xkno(2,kj)
                     yk = Xkno(2,kk)
!.....          Punkt ausserhalb Element:
                     ymin = MIN(yi,yj,yk)
                     IF ( Yp8>=ymin ) THEN
                        ymax = MAX(yi,yj,yk)
                        IF ( Yp8<=ymax ) THEN
!.....              Bis jetzt liegt Punkt innerhalb des das Element
!                   umschlieszenden Rechtecks XMIN/XMAX, YMIN/YMAX
!                   Pruefen, ob Punkt wirklich innerhalb DREIECK-Element
!                   liegt: BERECHNUNG DER TEILFLAECHEN (ohne 0.5)
                           f = xi*(yj-Yp8) + xj*(Yp8-yi) + Xp8*(yi-yj)
                           IF ( f>=0.0_rkind) THEN
                              f = xj*(yk-Yp8) + xk*(Yp8-yj) + Xp8*(yj-yk)
                              IF ( f>=0.0_rkind ) THEN
                                 f = xk*(yi-Yp8) + xi*(Yp8-yk) + Xp8*(yk-yi)
                                 IF ( f>=0.0_rkind ) THEN
                                    Ele = i
                                    ielem = Ele
                                    RETURN
                                 endif
                              endif
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
!     SCHLEIFE UEBER ALLE ELEMENTE wird hier folgendermassen hochgezaehlt:
!     beginnend bei IEALT, im Wechsel nach vorn und rueckwaerts suchend
         IF ( k<M .AND. if1==0 ) THEN
            if0 = 0
            IF ( l>1 ) if1 = 1
            k = k + 1
            i = k
            IF ( ijk<=M ) CYCLE
         endif
         CONTINUE
         EXIT
      ENDDO

      IF ( l>1 .AND. if0==0 ) THEN
         if1 = 0
         IF ( k<M ) if0 = 1
         l = l - 1
         i = l
      endif

      IF ( ijk<=M ) GOTO 100

      Ele = 0 
 
      END SUBROUTINE FIND_ELE_WIND
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function dintspec_y(ip,acloc,y)
      use datapool, only : msc,mdc,ddir,ds_incr,spsig,rkind, ZERO
      implicit none

      integer, intent(in)        :: ip
      real(rkind), intent(in)    :: y(msc), acloc(msc,mdc)

      integer             :: is, id
      real(rkind)         :: dintspec_y, tmp(msc)

      dintspec_y = ZERO
!     maxvalue   = maxval(ac2(:,:,ip))
!      if (maxvalue .lt. small) return 

      !acloc(:,:) = ac2(:,:,ip) !/ maxvalue

      do id = 1, mdc
        tmp(:) = acloc(:,id) * spsig * y
        do is = 2, msc
          dintspec_y = dintspec_y + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir 
        end do
      end do

      !dintspec_y = dintspec_y * maxvalue
          
      return
      end
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function dintspec(ip,acloc)
      use datapool,  only : msc,mdc,ddir,ds_incr,spsig,rkind, ZERO

      implicit none
      integer, intent(in)        :: ip
      real(rkind), intent(in)    :: acloc(msc,mdc)

      integer             :: is, id
      real(rkind)         :: dintspec, tmp(msc)

      dintspec = ZERO
      !maxvalue   = maxval(ac2(:,:,ip))
      !if (maxvalue .lt. small) return 

      !acloc(:,:) = ac2(:,:,ip) / maxvalue

      do id = 1, mdc
        tmp(:) = acloc(:,id) * spsig
        do is = 2, msc
          dintspec = dintspec+0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
        end do
      end do
      !dintspec = dintspec * maxvalue

      return
      end function
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTERDIR (Y1, Y2, DX, DIFFDX, YINTER)
      use datapool,  only : rkind, zero
      IMPLICIT NONE
      REAL(rkind), INTENT(IN)   :: Y1, Y2, DX, DIFFDX
      REAL(rkind), INTENT(OUT)  :: YINTER

      IF ((Y1 > ZERO .AND. Y1 < 90._rkind) .AND. (Y2 < 360._rkind .AND. Y2 > 270._rkind)) THEN
        YINTER=(Y1+360._rkind)+(Y2-(Y1+360._rkind))/DX*DIFFDX
        IF (YINTER > 360._rkind) YINTER = YINTER - 360._rkind
      ELSE IF  ((Y2 > ZERO .AND. Y2 < 90._rkind) .AND. (Y1 < 360._rkind .AND. Y1 > 270._rkind)) THEN
        YINTER = Y1+((Y2+360._rkind)-Y1)/DX*DIFFDX
        IF (YINTER > 360._rkind) YINTER = YINTER - 360._rkind
      ELSE
        YINTER = Y1+(Y2-Y1)/DX*DIFFDX
      END IF

      IF (ABS(YINTER) .LT. TINY(1.)) YINTER = Y1

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTERLIN (NX1, NX2, X1, X2, Y1, Y2)
      USE DATAPOOL, ONLY : THR, RKIND
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX1, NX2

      REAL(rkind), INTENT(IN)    :: X1(NX1), Y1(NX1)
      REAL(rkind), INTENT(IN)    :: X2(NX2)
      REAL(rkind), INTENT(OUT)   :: Y2(NX2)

      INTEGER             :: I, J
      REAL(rkind)                :: DX1(NX1-1)

      DO I = 1, NX1 - 1
        DX1(I) = X1(I+1) - X1(I)
      END DO

      DO I = 1, NX2
        DO J = 1, NX1 - 1
          IF (ABS(X2(I) - X1(J)) .LT. THR) THEN
            Y2(I) = Y1(J)
          ELSE IF (X2(I) .GT. X1(J) .AND. X2(I) .LT. X1(J+1)) THEN
            Y2(I) = Y1(J) + (Y1(J+1)-Y1(J))/DX1(J)*(X2(I)-X1(J))
          END IF
        END DO
      END DO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WWM_ABORT(string)
#ifdef MPI_PARALL_GRID
      USE DATAPOOL, only : parallel_abort, DBG
#else
      USE DATAPOOL, only : DBG
#endif

      IMPLICIT NONE
      character(*), intent(in) :: string

      WRITE(DBG%FHNDL, *) TRIM(string)
      FLUSH(DBG%FHNDL)

#ifdef MPI_PARALL_GRID
      CALL PARALLEL_ABORT(TRIM(string))
#else
      Print *, 'We have to abort. Reason:'
      Print *, TRIM(string)
      STOP 'WWM_ABORT'
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE TEST_FILE_EXIST_DIE(string1, string2)
      CHARACTER(LEN=*) :: string1
      CHARACTER(LEN=*) :: string2
      CHARACTER(LEN=512) :: ErrMsg
      LOGICAL :: LFLIVE
!      Print *, 'string1=', TRIM(string1)
!      Print *, 'string2=', TRIM(string2)
      INQUIRE( FILE = TRIM(string2), EXIST = LFLIVE )
      IF ( .NOT. LFLIVE ) THEN
        WRITE(ErrMsg,10) TRIM(string1), TRIM(string2)
  10    FORMAT(a, ' ', a)
        CALL WWM_ABORT(TRIM(ErrMsg))
      END IF
      END SUBROUTINE
!**********************************************************************
!* Simple Gauss method for solving systems                            *
!* We simply solve the equation A X = Y                               *
!**********************************************************************
      SUBROUTINE GAUSS_SOLVER(N, EMAT, X, Y)
      USE DATAPOOL, ONLY: RKIND, ONE
      IMPLICIT NONE
      integer, intent(in) :: N
      real(rkind), intent(inout) :: EMAT(N,N)
      real(rkind), intent(out) :: X(N)
      real(rkind), intent(in) :: Y(N)
      !
      real(rkind) :: EMATW(N,N), YW(N)
      integer U(N), Jsel, I, J, I2
      real(rkind) :: alpha, ThePivot, eVal, MaxVal
      U=0
      EMATW=EMAT
      YW=Y
      DO I=1,N
        MaxVal=0
        Jsel=-1
        DO J=1,N
          IF (U(J) .eq. 0) THEN
            eVal=EMATW(J,I)
            IF (abs(eVal) .gt. MaxVal) THEN
              MaxVal=abs(eVal)
              Jsel=J
            END IF
          END IF
        END DO
        IF (Jsel == -1) THEN
          CALL WWM_ABORT('Error in GAUSS_SOLVER, singular matrix')
        END IF
        ThePivot=ONE/EMATW(Jsel,I)
        U(Jsel)=I
        DO I2=I,N
          EMATW(Jsel,I2)=EMATW(Jsel,I2)*ThePivot
        END DO
        YW(Jsel)=YW(Jsel)*ThePivot
        DO J=1,N
          IF (J .ne. Jsel) THEN
            alpha=EMATW(J,I)
            DO I2=I,N
              EMATW(J,I2)=EMATW(J,I2) - alpha*EMATW(Jsel,I2)
            END DO
            YW(J)=YW(J) - alpha*YW(Jsel)
          END IF
        END DO
      END DO
      DO Jsel=1,N
        I=U(Jsel)
        X(I)=YW(Jsel)
      END DO
      END SUBROUTINE
!**********************************************************************
!* from http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm     *
!**********************************************************************
      subroutine SOLVE_TRIDIAG(a,b,c,d,x,n)
      use datapool, only : rkind
      implicit none
!        a - sub-diagonal (means it is the diagonal below the main diagonal)
!        b - the main diagonal
!        c - sup-diagonal (means it is the diagonal above the main diagonal)
!        d - right part
!        x - the answer
!        n - number of equations
      integer,intent(in) :: n
      real(rkind), intent(in) :: a(n), b(n), c(n), d(n)
      real(rkind), intent(out) :: x(n)
      real(rkind) :: cp(n),dp(n)
      real(rkind) :: m
      integer i
 
! initialize c-prime and d-prime
      cp(1) = c(1)/b(1)
      dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
      do i = 2,n
        m = b(i)-cp(i-1)*a(i)
        cp(i) = c(i)/m
        dp(i) = (d(i)-dp(i-1)*a(i))/m
      enddo
! initialize x
      x(n) = dp(n)
! solve for x from the vectors c-prime and d-prime
      do i = n-1, 1, -1
        x(i) = dp(i)-cp(i)*x(i+1)
      end do
      end subroutine SOLVE_TRIDIAG
!**********************************************************************
!*                                                                    *
!**********************************************************************
            SUBROUTINE GAUS1D( M, AMAT, R, AM1, A1M )
               USE DATAPOOL, ONLY: RKIND
               IMPLICIT NONE

               INTEGER    :: M, I, K
               REAL(rkind)       :: AMAT(3,M), R(M), A1M, AM1, FAK
               REAL(rkind)       :: SPALTE(M), ZEILE(M)

! AMAT: TRIDIAGONAL: 1. LEFT DIAGONALE; 2. MAIN DIAGONAL; 3. RIGHT DIAGONAL
! DIMENSION OF MATRIX
! AM1 LOWER  LEFT ELEMENT
! A1M UPPER RIGHT ELEMENT
! M ~ EQUATIONS
! R ~ RIGHT HAND SIDE

               DO K = 1, M
                  ZEILE(K) = 0.0
                  SPALTE(K)= 0.0
               END DO
               SPALTE(1) = A1M
               SPALTE(M-1) = AMAT(3,M-1)
               SPALTE(M) = AMAT(2,M)
               ZEILE(1) = AM1
! R
               DO I = 2, M-1
                  FAK = AMAT(1,I)/AMAT(2,I-1)
                  R(I)= R(I) - R(I-1)*FAK
                  AMAT(2,I) = AMAT(2,I) - AMAT(3,I-1)*FAK
                  SPALTE(I) = SPALTE(I) - SPALTE(I-1)*FAK
               END DO
!               I = M
               ZEILE(M-1) = AMAT(1,M)

               DO K = 2, M-1
                  FAK = ZEILE(K-1) / AMAT(2,K-1)
                  ZEILE(K) = ZEILE(K) - AMAT(3,K-1)*FAK
                  R(M) = R(M) - R(K-1)*FAK
                  SPALTE(M) = SPALTE(M) - SPALTE(K-1)*FAK
               END DO
               K = M
               FAK = ZEILE(K-1) / AMAT(2,K-1)
               R(M) = R(M) - R(K-1)*FAK
               SPALTE(M) = SPALTE(M) - SPALTE(K-1)*FAK
! R
               I = M - 1
               FAK = SPALTE(I) / SPALTE(M)
               R(I) = R(I) - R(M)*FAK
               AMAT(3,I) = 0.0
               DO I = M-2, 1, -1
                  FAK = SPALTE(I) / SPALTE(M)
                  R(I) = R(I) - R(M)*FAK
                  AMAT(3,I) = AMAT(3,I) - R(M)*FAK
                  FAK = AMAT(3,I) / AMAT(2,I+1)
                  R(I) = R(I) - R(I+1)*FAK
               END DO

               DO I = 1, M-1
                  R(I) = R(I) / AMAT(2,I)
               END DO

               R(M) = R(M) / SPALTE(M)
            END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************

      SUBROUTINE SSORT2 (X, Y, Z, N, KFLAG)
      USE DATAPOOL, only : rkind
!TS: Modified second array z(*) to carry along

!***BEGIN PROLOGUE  SSORT
!***PURPOSE  Sort an array and optionally make the same interchanges in
!            an auxiliary array.  The array may be sorted in increasing
!            or decreasing order.  A slightly modified QUICKSORT
!            algorithm is used.
!***LIBRARY   SLATEC
!***CATEGORY  N6A2B
!***TYPE      SINGLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
!***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
!***AUTHOR  Jones, R. E., (SNLA)
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!
!   SSORT sorts array X and optionally makes the same interchanges in
!   array Y.  The array X may be sorted in increasing order or
!   decreasing order.  A slightly modified quicksort algorithm is used.
!
!   Description of Parameters
!      X - array of values to be sorted   (usually abscissas)
!      Y - array to be (optionally) carried along
!      N - number of values in array X to be sorted
!      KFLAG - control parameter
!            =  2  means sort X in increasing order and carry Y along.
!            =  1  means sort X in increasing order (ignoring Y)
!            = -1  means sort X in decreasing order (ignoring Y)
!            = -2  means sort X in decreasing order and carry Y along.
!
!***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
!                 for sorting with minimal storage, Communications of
!                 the ACM, 12, 3 (1969), pp. 185-187.
!***REVISION HISTORY  (YYMMDD)
!   761101  DATE WRITTEN
!   761118  Modified to use the Singleton quicksort algorithm.  (JAW)
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891009  Removed unreferenced statement labels.  (WRB)
!   891024  Changed category.  (WRB)
!   891024  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   901012  Declared all variables; changed X,Y to SX,SY. (M. McClain)
!   920501  Reformatted the REFERENCES section.  (DWL, WRB)
!   920519  Clarified error messages.  (DWL)
!   920801  Declarations section rebuilt and code restructured to use
!           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
!***END PROLOGUE  SSORT
!     .. Scalar Arguments ..
      INTEGER KFLAG, N
!     .. Array Arguments ..
      REAL(rkind) X(*), Y(*), Z(*)
!     .. Local Scalars ..
      REAL(rkind) R, T, TT, TTY, TY, TZ, TTZ
      INTEGER I, IJ, J, K, KK, L, M, NN
!     .. Local Arrays ..
      INTEGER IL(21), IU(21)
!     .. External Subroutines ..
!     None
!     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
!***FIRST EXECUTABLE STATEMENT  SSORT
      NN = N
      IF (NN .LT. 1) THEN
        WRITE (*,*) 'The number of values to be sorted is not positive.'
         RETURN
      ENDIF
!
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
         WRITE (*,*) 'The sort control parameter, K, is not 2, 1, -1, or -2.'
         RETURN
      ENDIF
!
!     Alter array X to get decreasing order if needed
!
      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            X(I) = -X(I)
   10    CONTINUE
      ENDIF
!
      IF (KK .EQ. 2) GO TO 100
!
!     Sort X only
!
      M = 1
      I = 1
      J = NN
      R = 0.375E0
!
   20 IF (I .EQ. J) GO TO 60
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
!
   30 K = I
!
!     Select a central element of the array and save it in location T
!
      IJ = I + INT((J-I)*R)
      T = X(IJ)
!
!     If first element of array is greater than T, interchange with T
!
      IF (X(I) .GT. T) THEN
         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)
      ENDIF
      L = J
!
!     If last element of array is less than than T, interchange with T
!
      IF (X(J) .LT. T) THEN
         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)
!
!        If first element of array is greater than T, interchange with T
!
         IF (X(I) .GT. T) THEN
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
         ENDIF
      ENDIF
!
!     Find an element in the second half of the array which is smaller
!     than T
!
   40 L = L-1
      IF (X(L) .GT. T) GO TO 40
!
!     Find an element in the first half of the array which is greater
!     than T
!
   50 K = K+1
      IF (X(K) .LT. T) GO TO 50
!
!     Interchange these elements
!
      IF (K .LE. L) THEN
         TT = X(L)
         X(L) = X(K)
         X(K) = TT
         GO TO 40
      ENDIF
!
!     Save upper and lower subscripts of the array yet to be sorted
!
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 70
!
!     Begin again on another portion of the unsorted array
!
   60 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
!
   70 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1
!
   80 I = I+1
      IF (I .EQ. J) GO TO 60
      T = X(I+1)
      IF (X(I) .LE. T) GO TO 80
      K = I
!
   90 X(K+1) = X(K)
      K = K-1
      IF (T .LT. X(K)) GO TO 90
      X(K+1) = T
      GO TO 80
!
!     Sort X and carry Y along
!
  100 M = 1
      I = 1
      J = NN
      R = 0.375E0
!
  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
!
  120 K = I
!
!     Select a central element of the array and save it in location T
!
      IJ = I + INT((J-I)*R)
      T = X(IJ)
      TY = Y(IJ)
      TZ = Z(IJ)
!
!     If first element of array is greater than T, interchange with T
!
      IF (X(I) .GT. T) THEN
         
         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)
         
         Y(IJ) = Y(I)
         Y(I) = TY
         TY = Y(IJ)
         
         Z(IJ) = Z(I)
         Z(I) = TZ
         TZ = Z(IJ)
      ENDIF
      L = J
!
!     If last element of array is less than T, interchange with T
!
      IF (X(J) .LT. T) THEN
      
         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)
      
         Y(IJ) = Y(J)
         Y(J) = TY
         TY = Y(IJ)
         
         Z(IJ) = Z(J)
         Z(J) = TZ
         TZ = Z(IJ)
!
!        If first element of array is greater than T, interchange with T
!
         IF (X(I) .GT. T) THEN
            
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
            
            Y(IJ) = Y(I)
            Y(I) = TY
            TY = Y(IJ)
            
            Z(IJ) = Z(I)
            Z(I) = TZ
            TZ = Z(IJ)
            
         ENDIF
      ENDIF
!
!     Find an element in the second half of the array which is smaller
!     than T
!
  130 L = L-1
      IF (X(L) .GT. T) GO TO 130
!
!     Find an element in the first half of the array which is greater
!     than T
!
  140 K = K+1
      IF (X(K) .LT. T) GO TO 140
!
!     Interchange these elements
!
      IF (K .LE. L) THEN
      
         TT = X(L)
         X(L) = X(K)
         X(K) = TT
      
         TTY = Y(L)
         Y(L) = Y(K)
         Y(K) = TTY
      
         TTZ = Z(L)
         Z(L) = Z(K)
         Z(K) = TTZ
      
         GO TO 130
      ENDIF
!
!     Save upper and lower subscripts of the array yet to be sorted
!
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160
!
!     Begin again on another portion of the unsorted array
!
  150 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
!
  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1
!
  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = X(I+1)
      TY = Y(I+1)
      TZ = Z(I+1)
      IF (X(I) .LE. T) GO TO 170
      K = I
!
  180 X(K+1) = X(K)
      Y(K+1) = Y(K)
      Z(K+1) = Z(K)
      K = K-1
      IF (T .LT. X(K)) GO TO 180
      X(K+1) = T
      Y(K+1) = TY
      Z(K+1) = TZ
      GO TO 170
!
!     Clean up
!
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            X(I) = -X(I)
  200    CONTINUE
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                     *
!**********************************************************************
      FUNCTION POSITION_BEFORE_POINT(FILEIN)
      IMPLICIT NONE
      character(len=*), intent(in) :: FILEIN
      character(len=1), parameter :: ePoint = '.'
      integer POSITION_BEFORE_POINT
      integer len, LPOS, I
      len=LEN_TRIM(FILEIN)
      LPOS=-1
      DO I=1,len
        IF (FILEIN(I:I).eq.ePoint) THEN
          LPOS=I-1
        END IF
      ENDDO
      IF (LPOS.eq.-1) THEN
        LPOS=len
      ENDIF
      POSITION_BEFORE_POINT=LPOS
      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GETSTRING (ThePow, TheNb, eStr)
      IMPLICIT NONE
      integer, intent(in) :: TheNb, ThePow
      character(len=ThePow), intent(out) :: eStr
      character(len=ThePow) :: eStrZero
      character(len=40) :: eStrNb
      INTEGER :: ePower, iPower, WeFind, I, nb, IFIRST, h
      IF (TheNb.lt.0) THEN
        CALL WWM_ABORT('GETSTRING - call with negative value')
      END IF
      ePower=1
      WeFind=0
      DO I=1,ThePow
        ePower=ePower*10
        IF ((WeFind.eq.0).and.(TheNb.lt.ePower)) THEN
          WeFind=1
          iPower=I
        END IF
      END DO
      IF (WeFind.eq.0) THEN
        CALL WWM_ABORT('GETSTRING - call with non-representable number')
      END IF
      nb=ThePow - iPower
      DO I=1,ThePow
        eStrZero(I:I)=' '
      END DO
      DO I=1,nb
        eStrZero(I:I)='0'
      END DO
      WRITE(eStrNb,*) TheNb
      IFIRST=-1
      DO I=1,40
        IF ((eStrNb(I:I).ne.' ').and.(IFIRST.eq.-1)) THEN
          IFIRST=I
        ENDIF
      END DO
      nb=40+1-IFIRST
      DO I=1,nb
        h=I+IFIRST-1
        eStrNb(I:I)=eStrNb(h:h)
      END DO
      DO I=nb+1,40
        eStrNb(I:I)=' '
      END DO
      WRITE(eStr,40) TRIM(eStrZero),TRIM(eStrNb)
  40  FORMAT (a,a)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
     SUBROUTINE SPECSTAT(TIME)

         USE DATAPOOL
         IMPLICIT NONE

         REAL, INTENT(IN)    :: TIME

         INTEGER             :: IS, ID, IP
         REAL(RKIND)         :: ETOT

! Specstat 2D definiton 1 based on the total energy ... may be to less specific for certain spatial points e.g. coastal points in shadowing areas ....

        STAT2D = 0.
        DO IP = 1, MNP
           ETOT = SUM(AC2(:,:,IP))
           DO IS = 1, MSC
             DO ID = 1, MDC
               IF (ETOT .GT. THR) THEN
                 IF (AC2(IS,ID,IP)/ETOT .LT. 1.E-5) THEN
                   STAT2D(IS,ID) = STAT2D(IS,ID) + 1./REAL(MNP)*100.
!                   WRITE(*,'(3I10,4E15.4)') IP, IS, ID, STAT2D(IS,ID), AC2(IS,ID,IP)/ETOT, AC2(IS,ID,IP), ETOT
                 END IF
               ELSE 
                 STAT2D(IS,ID) = STAT2D(IS,ID) + 1./REAL(MNP)*100.
               END IF
             END DO
           END DO
        END DO

        DO IS = 1, MSC
          DO ID = 1, MDC
            IF (STAT2D(IS,ID) .GT. 99.9) STAT2D(IS,ID) = 0.
          END DO
        END DO

!         CALL SSORTORG(STAT1D,SIGTMP,MSC,1)
!         DO IS = 1, MSC
!            DO ID = 1, MDC
!              WRITE(*,'(2I10,F15.4)') IS, ID, STAT2D(IS,ID)
!            END DO
!             WRITE(*,*) IS, STAT1D(IS)
!          END DO

!         CALL DISPLAY_GRAPH(SPSIG,SUM1D,MSC,MINVAL(SPSIG),MAXVAL(SPSIG),MINVAL(SUM1D),MAXVAL(SUM1D),'','','')

       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE INTER_STRUCT_DATA(NDX,NDY,DX,DY,OFFSET_X,OFFSET_Y,MAT,VAL)
       USE DATAPOOL
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: NDX, NDY
       REAL(rkind), INTENT(IN)    :: MAT(NDX,NDY)
       REAL(rkind), INTENT(IN)    :: DX, DY, OFFSET_X, OFFSET_Y
       REAL(rkind), INTENT(OUT)   :: VAL(MNP_WIND)
       INTEGER             :: IP, J_INT, I_INT
       REAL(rkind)                :: WX1, WX2, WX3, WX4, HX1, HX2
       REAL(rkind)                :: DELTA_X, DELTA_Y, LEN_X, LEN_Y
       DO IP = 1, MNP_WIND
         LEN_X = XP_WIND(IP)-OFFSET_X
         LEN_Y = YP_WIND(IP)-OFFSET_Y
         I_INT = INT( LEN_X/DX ) + 1
         J_INT = INT( LEN_Y/DY ) + 1
         DELTA_X = LEN_X - (I_INT - 1) * DX ! Abstand X u. Y
         DELTA_Y = LEN_Y - (J_INT - 1) * DY !
         WX1     = MAT(  I_INT   , J_INT  ) ! Unten Links
         WX2     = MAT(  I_INT   , J_INT+1) ! Oben  Links
         WX3     = MAT(  I_INT+1,  J_INT+1) ! Oben  Rechts
         WX4     = MAT(  I_INT+1,  J_INT  ) ! Unten Rechts
         IF (IWINDFORMAT == 2) THEN
           HX1 = WX1 + (WX2-WX1)/DX * DELTA_X
           HX2 = WX4 + (WX3-WX4)/DX * DELTA_X
         ELSE IF (IWINDFORMAT == 3) THEN
           HX1 = WX1 + (WX4-WX1)/DX * DELTA_X
           HX2 = WX2 + (WX3-WX2)/DX * DELTA_X
         ELSE
           CALL WWM_ABORT('Write your HX1/HX2 code here')
         END IF
         VAL(IP) = HX1 + (HX2-HX1)/DY * DELTA_Y
       END DO
       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CreateAngleMatrix(eta_rho, xi_rho, ANG_rho, LON_rho, LAT_rho)
      implicit none
      integer, intent(in) :: eta_rho, xi_rho
      REAL*8, DIMENSION(eta_rho, xi_rho), intent(in) :: LON_rho, LAT_rho
      REAL*8, DIMENSION(eta_rho, xi_rho), intent(out) :: ANG_rho
      !
      integer eta_u, xi_u, iEta, iXi
      real*8, allocatable :: LONrad_u(:,:)
      real*8, allocatable :: LATrad_u(:,:)
      real*8, allocatable :: azim(:,:)
      real*8 :: eAzim, fAzim, dlam, eFact1, eFact2
      real*8 :: signAzim, signDlam, ThePi, DegTwoRad
      real*8 :: eLon, eLat, phi1, phi2, xlam1, xlam2
      real*8 :: TPSI2, cta12
      integer istat
      eta_u=eta_rho
      xi_u=xi_rho-1
      allocate(LONrad_u(eta_u,xi_u), LATrad_u(eta_u,xi_u), azim(eta_u,xi_u-1), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 33')
      ThePi=3.141592653589792
      DegTwoRad=ThePi/180
      DO iEta=1,eta_u
        DO iXi=1,xi_u
          eLon=(LON_rho(iEta,iXi)+LON_rho(iEta,iXi+1))*0.5
          eLat=(LAT_rho(iEta,iXi)+LAT_rho(iEta,iXi+1))*0.5
          LONrad_u(iEta,iXi)=eLon*DegTwoRad
          LATrad_u(iEta,iXi)=eLat*DegTwoRad
        END DO
      END DO
      DO iEta=1,eta_u
        DO iXi=1,xi_u-1
          phi1=LATrad_u(iEta,iXi)
          xlam1=LONrad_u(iEta,iXi)
          phi2=LATrad_u(iEta,iXi+1)
          xlam2=LONrad_u(iEta,iXi+1)
          TPSI2=TAN(phi2)
          dlam=xlam2-xlam1
          CALL TwoPiNormalization(dlam)
          cta12=(cos(phi1)*TPSI2 - sin(phi1)*cos(dlam))/sin(dlam)
          eAzim=ATAN(1./cta12)
          CALL MySign(eAzim, signAzim)
          CALL MySign(dlam, signDlam)
          IF (signDlam.ne.signAzim) THEN
            eFact2=1
          ELSE
            eFact2=0
          END IF
          eFact1=-signAzim
          fAzim=eAzim+ThePi*eFact1*eFact2
          azim(iEta,iXi)=fAzim
        END DO
      END DO
      DO iEta=1,eta_u
        DO iXi=2,xi_u
          ANG_rho(iEta,iXi)=ThePi*0.5 - azim(iEta,iXi-1)
        END DO
      END DO
      DO iEta=1,eta_u
        ANG_rho(iEta,1)=ANG_rho(iEta,2)
        ANG_rho(iEta,xi_rho)=ANG_rho(iEta,xi_u)
      END DO
      deallocate(LONrad_u, LATrad_u, azim)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CreateAngleMatrix_v(eta_rho,xi_rho,ANG_rho,LON_rho,LAT_rho)
      USE DATAPOOL, only : rkind, istat
      implicit none
      integer, intent(in) :: eta_rho, xi_rho
      REAL(rkind), DIMENSION(eta_rho, xi_rho), intent(in) :: LON_rho, LAT_rho
      REAL(rkind), DIMENSION(eta_rho, xi_rho), intent(out) :: ANG_rho
      !
      integer eta_v, xi_v, iEta, iXi
      real(rkind), allocatable :: LONrad_v(:,:)
      real(rkind), allocatable :: LATrad_v(:,:)
      real(rkind), allocatable :: azim(:,:)
      real(rkind) :: eAzim, fAzim, dlam, eFact1, eFact2
      real(rkind) :: signAzim, signDlam, ThePi, DegTwoRad
      real(rkind) :: eLon, eLat, phi1, phi2, xlam1, xlam2
      real(rkind) :: TPSI2, cta12
      eta_v=eta_rho-1
      xi_v=xi_rho
      allocate(LONrad_v(eta_v,xi_v), LATrad_v(eta_v,xi_v), azim(eta_v-1,xi_v), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 34')

      ThePi=3.141592653589792
      DegTwoRad=ThePi/180
      DO iEta=1,eta_v
        DO iXi=1,xi_v
          eLon=(LON_rho(iEta,iXi)+LON_rho(iEta+1,iXi))*0.5
          eLat=(LAT_rho(iEta,iXi)+LAT_rho(iEta+1,iXi))*0.5
          LONrad_v(iEta,iXi)=eLon*DegTwoRad
          LATrad_v(iEta,iXi)=eLat*DegTwoRad
        END DO
      END DO
      DO iEta=1,eta_v-1
        DO iXi=1,xi_v
          phi1=LATrad_v(iEta,iXi)
          xlam1=LONrad_v(iEta,iXi)
          phi2=LATrad_v(iEta+1,iXi)
          xlam2=LONrad_v(iEta+1,iXi)
          TPSI2=TAN(phi2)
          dlam=xlam2-xlam1
          CALL TwoPiNormalization(dlam)
          cta12=(cos(phi1)*TPSI2 - sin(phi1)*cos(dlam))/sin(dlam)
          eAzim=ATAN(1./cta12)
          CALL MySign(eAzim, signAzim)
          CALL MySign(dlam, signDlam)
          IF (signDlam.ne.signAzim) THEN
            eFact2=1
          ELSE
            eFact2=0
          END IF
          eFact1=-signAzim
          fAzim=eAzim+ThePi*eFact1*eFact2
          azim(iEta,iXi)=fAzim
        END DO
      END DO
      DO iEta=2,eta_v
        DO iXi=1,xi_v
          ANG_rho(iEta,iXi)=ThePi*0.5 - azim(iEta-1,iXi)
        END DO
      END DO
      DO iXi=1,xi_v
        ANG_rho(1,iXi)=ANG_rho(2,iXi)
        ANG_rho(eta_rho,iXi)=ANG_rho(eta_v,iXi)
      END DO
      deallocate(LONrad_v)
      deallocate(LATrad_v)
      deallocate(azim)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE TwoPiNormalization(TheAng)
      USE DATAPOOL, only : rkind
      implicit none
      REAL(rkind), intent(inout) :: TheAng
      !
      REAL(rkind) :: ThePi
      ThePi=3.141592653589792
      IF (TheAng < -ThePi) THEN
        TheAng=TheAng + 2*ThePi
      ENDIF
      IF (TheAng > ThePi) THEN
        TheAng=TheAng - 2*ThePi
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MySign(TheVal, TheSign)
      USE DATAPOOL, only : rkind
      implicit none
      REAL(rkind), intent(in) :: TheVal
      REAL(rkind), intent(out) :: TheSign
      IF (TheVal > 0) THEN
        TheSign=1
      ELSEIF (TheVal < 0) THEN
        TheSign=-1
      ELSE
        TheSign=0
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
