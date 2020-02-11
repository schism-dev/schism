#include "wwm_functions.h"
!#define DEBUG
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SPHERICAL_COORDINATE_DISTANCE(PT1_LON, PT2_LON, PT1_LAT, PT2_LAT, DIST)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(RKIND), intent(in)  :: PT1_LON, PT2_LON, PT1_LAT, PT2_LAT
      REAL(RKIND), intent(out) :: DIST
      REAL(RKIND) :: PT1_LON_RAD, PT2_LON_RAD, PT1_LAT_RAD, PT2_LAT_RAD
      REAL(RKIND) :: COEFF
      REAL(RKIND) :: X1, Y1, Z1, X2, Y2, Z2
      REAL(RKIND) :: SCALPROD

      ! Point coordinates in a Cartesian grid
      COEFF = PI/180.0_rkind

      ! Computing the distance
      PT1_LON_RAD = PT1_LON*COEFF
      PT2_LON_RAD = PT2_LON*COEFF
      PT1_LAT_RAD = PT1_LAT*COEFF
      PT2_LAT_RAD = PT2_LAT*COEFF
      X1 = COS(PT1_LON_RAD)*COS(PT1_LAT_RAD)
      Y1 = SIN(PT1_LON_RAD)*COS(PT1_LAT_RAD)
      Z1 = SIN(PT1_LAT_RAD)
      X2 = COS(PT2_LON_RAD)*COS(PT2_LAT_RAD)
      Y2 = SIN(PT2_LON_RAD)*COS(PT2_LAT_RAD)
      Z2 = SIN(PT2_LAT_RAD)
      SCALPROD = X1*X2+Y1*Y2+Z1*Z2;
      if (SCALPROD .GE. 1) THEN
        DIST = 0;
      ELSE
        DIST = ACOS(SCALPROD)
      END IF
      END SUBROUTINE
!********************************************************************************
!* See formula for area of spherical triangle in                                *
!* http://math.stackexchange.com/questions/9819/area-of-a-spherical-triangle    *
!*                                                                              *
!********************************************************************************
      SUBROUTINE SPHERICAL_COORDINATE_AREA(LON1, LON2, LON3, LAT1, LAT2, LAT3, AREA)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), intent(in) :: LON1, LON2, LON3, LAT1, LAT2, LAT3
      REAL(rkind), intent(out) :: AREA
      REAL(rkind) :: DistA, DistB, DistC, DistS
      REAL(rkind) :: eTan1, eTan2, eTan3, eTan4
      REAL(rkind) :: eProd, sqrtProd
      CALL SPHERICAL_COORDINATE_DISTANCE(LON1, LON2, LAT1, LAT2, DistA)
      CALL SPHERICAL_COORDINATE_DISTANCE(LON1, LON3, LAT1, LAT3, DistB)
      CALL SPHERICAL_COORDINATE_DISTANCE(LON2, LON3, LAT2, LAT3, DistC)
      DistS=(DistA + DistB + DistC)/2.0_rkind
      eTan1=tan(DistS/2.0_rkind)
      eTan2=tan((DistS - DistA)/2.0_rkind)
      eTan3=tan((DistS - DistB)/2.0_rkind)
      eTan4=tan((DistS - DistC)/2.0_rkind)
      eProd=eTan1*eTan2*eTan3*eTan4
      sqrtProd=SQRT(eProd)
      AREA=4.0_rkind*ATAN(sqrtProd)
      END SUBROUTINE


      SUBROUTINE CORRECT_SINGLE_DXP(DXP)
      USE DATAPOOL
      IMPLICIT NONE
      real(rkind), intent(inout) :: DXP
      IF (DXP .le. -180_rkind) THEN
        DXP=DXP + 360
      END IF
      IF (DXP .ge. 180) THEN
        DXP=DXP - 360
      END IF
      END SUBROUTINE

      SUBROUTINE INIT_SPATIAL_GRID
         USE DATAPOOL
         IMPLICIT NONE
         REAL(rkind) :: TL1, TL2, TL3
         REAL(rkind) :: TMPTLMIN
         REAL(rkind) :: TMPTLMAX 
         REAL(rkind) :: DBLTMP, DXP1, DXP2, DXP3, DYP1, DYP2, DYP3
         REAL(rkind) :: P1_XLOC, P2_XLOC, P3_XLOC
         REAL(rkind) :: PROV1, PROV2, PROV3
         REAL(rkind) :: AREA, AREA_RAD
         INTEGER           :: I1, I2, I3, TMPINE, NI(3)
         INTEGER           :: J1, J2, J3
         INTEGER           :: IP, IE, IEWRONG, IEWRONGSUM
         LOGICAL           :: LWRONG

#ifdef MPI_PARALL_GRID
         REAL(rkind)            :: AVETA_GL, TLMIN_GL, TLMAX_GL
         AVETA_GL = ZERO; TLMIN_GL = ZERO; TLMAX_GL = ZERO
#endif

         TMPTLMIN = 10.E10_rkind
         TMPTLMAX = ZERO
         TLMIN = VERYLARGE 
         TLMAX = VERYSMALL 
         AVETL = ZERO
         AVETA = ZERO

         LWRONG = .FALSE.

         SELECT CASE (DIMMODE)
!
!     *** One dimension mode
!
             CASE (1)

                IF (LVAR1D) THEN
                  DX1(0)     = XP(2)- XP(1)
                  DX1(1)     = DX1(0)
                  DX1(MNP)   = XP(MNP) - XP(MNP-1)
                  DX1(MNP+1) = DX1(MNP)
                  DX2(0)     = DX1(0)
                  DX2(MNP+1) = DX1(MNP)
                  DO IP = 2, MNP-1 ! Bandwith at gridpoints
                     DX1(IP) = (XP(IP)-XP(IP-1))/2. + (XP(IP+1)-XP(IP))/2.
                  END DO
                  DO IP = 2, MNP ! Stepwidth between gridpoints K and K-1
                     DX2(IP) = XP(IP) - XP(IP-1)
                  END DO
                  DX2(1) = DX1(0)
                END IF

                TL1 = 0.
                DO IE = 1, MNP-1
                   TL1 = TL1 + ABS(XP(IE+1) - XP(IE))
                END DO
                AVETL = TL1/MyREAL(MNP)
!
!     *** Two dimension mode
!
             CASE(2)

                IEWRONGSUM = 0
                DO IE = 1, MNE

                   I1 = INE(1,IE)
                   I2 = INE(2,IE)
                   I3 = INE(3,IE)
                   NI = INE(:,IE)

                   ! Special treatment for elements overpassing the dateline
                   P1_XLOC = XP(I1);
                   P2_XLOC = XP(I2);
                   P3_XLOC = XP(I3);
                   IF (LSPHE) THEN
                     ! Kevin Martins' fix to close the mesh at the
                     ! dateline: understand whether an element has a single node
                     ! isolated on the other side, and bring it
                     ! back. This approach does nothing with the element
                     ! that contains the pole, and apparently is stable
                     ! there.
                     IF (    (P1_XLOC - P2_XLOC).GT.180 .AND. & 
                           & (P3_XLOC - P2_XLOC).GT.180 ) THEN
                       ! In this case, P2 is 'isolated' to the East of the dateline; we "bring it back" towards the West
                       P2_XLOC = 180.0_rkind + ABS(-180.0_rkind - P2_XLOC)
                     ELSE IF (    (P2_XLOC - P1_XLOC).GT.180 .AND. &
                                & (P3_XLOC - P1_XLOC).GT.180  ) THEN
                       ! In this case, P1 is 'isolated' to the East of the dateline; we "bring it back" towards the West
                       P1_XLOC = 180.0_rkind + ABS(-180.0_rkind - P1_XLOC)
                     ELSE IF (    (P1_XLOC - P3_XLOC).GT.180 .AND. &
                                & (P2_XLOC - P3_XLOC).GT.180  ) THEN
                       ! In this case, P3 is 'isolated' to the East of the dateline; we "bring it back" towards the West
                       P3_XLOC = 180.0_rkind + ABS(-180.0_rkind - P3_XLOC)
                     ELSE IF (    (P1_XLOC - P2_XLOC).GT.180 .AND. &  
                                & (P1_XLOC - P3_XLOC).GT.180  ) THEN
                       ! In this case, P1 is 'isolated' to the West of the dateline; we "bring it back" towards the East
                       P1_XLOC = -180.0_rkind - ABS(180.0_rkind - P1_XLOC)
                     ELSE IF (    (P2_XLOC - P1_XLOC).GT.180 .AND. &  
                                & (P2_XLOC - P3_XLOC).GT.180  ) THEN
                       ! In this case, P2 is 'isolated' to the West of the dateline; we "bring it back" towards the East
                       P2_XLOC = -180.0_rkind - ABS(180.0_rkind - P2_XLOC)
                     ELSE IF (    (P3_XLOC - P1_XLOC).GT.180 .AND. &  
                                & (P3_XLOC - P2_XLOC).GT.180  ) THEN
                       ! In this case, P3 is 'isolated' to the West of the dateline; we "bring it back" towards the East
                       P3_XLOC = -180.0_rkind - ABS(180.0_rkind - P3_XLOC)
                     END IF
                   END IF

                   IF (IGRIDTYPE.ne.2) THEN
                     DXP1 = P2_XLOC - P1_XLOC
                     DYP1 = YP(I2) - YP(I1)
                     DXP2 = P3_XLOC - P2_XLOC
                     DYP2 = YP(I3) - YP(I2)
                     DXP3 = P1_XLOC - P3_XLOC
                     DYP3 = YP(I1) - YP(I3)

                     IF (APPLY_DXP_CORR) THEN
                       ! the option APPLY_DXP_CORR does something
                       ! similar to Kevin Martin's approach to close the
                       ! mesh at the dateline. If the longitudinal
                       ! distance DXP is >180deg or <-180deg it is
                       ! adjusted subtracting or adding 360deg. This
                       ! approach also acts on the element containing
                       ! the pole, correcting a single longitudinal
                       ! disnance of the three, making this element
                       ! unstable.
                       ! This option could be a substitute of Kevin
                       ! Martin's fix if a special treatment of the
                       ! element containing the pole is introduced.
                       ! This option could be removed if it is a
                       ! duplication of Kevin Martin's fix.
                       CALL CORRECT_SINGLE_DXP(DXP1)
                       CALL CORRECT_SINGLE_DXP(DXP2)
                       CALL CORRECT_SINGLE_DXP(DXP3)
                     END IF

                     IEN(1,IE) = - DYP2
                     IEN(2,IE) =   DXP2
                     IEN(3,IE) = - DYP3
                     IEN(4,IE) =   DXP3
                     IEN(5,IE) = - DYP1
                     IEN(6,IE) =   DXP1
                     DBLTMP = (DXP3*DYP1 - DYP3*DXP1)*ONEHALF
                     IF (LSPHE .and. USE_EXACT_FORMULA_SPHERICAL_AREA) THEN
                       CALL SPHERICAL_COORDINATE_AREA(XP(I1), XP(I2), XP(I3), YP(I1), YP(I2), YP(I3), AREA_RAD)
                       AREA=AREA_RAD*RADDEG*RADDEG
#ifdef DEBUG
                       WRITE(STAT%FHNDL,*) 'IE=', IE
                       WRITE(STAT%FHNDL,*) 'I123=', I1, I2, I3
                       WRITE(STAT%FHNDL,*) 'XP123=', XP(I1), XP(I2), XP(I3)
                       WRITE(STAT%FHNDL,*) 'YP123=', YP(I1), YP(I2), YP(I3)
                       WRITE(STAT%FHNDL,*) '  AREA=', AREA
                       WRITE(STAT%FHNDL,*) 'DBLTMP=', DBLTMP
                       FLUSH(STAT%FHNDL)
#endif
                     ELSE
                       AREA=DBLTMP
                     END IF
                     TRIA(IE) = AREA
                   END IF


                   IF (TRIA(IE) .LT. 0.0) THEN
                      TMPINE = INE(2,IE)
                      INE(2,IE) = INE(3,IE)
                      INE(3,IE) = TMPINE
                      IF (ANY(INE(:,IE) .GT. MNP)) CALL WWM_ABORT('WRITE ELEMENT CONNECTION TALBE HAS NODENUMBERS GT MNP')
                      I2 = INE(2,IE)
                      I3 = INE(3,IE)
                      TRIA(IE) = -TRIA(IE)
                      PROV1=IEN(6,IE) ! DXP1
                      PROV2=IEN(2,IE) ! DXP2
                      PROV3=IEN(4,IE) ! DXP3
                      IEN(6,IE)=-PROV3
                      IEN(2,IE)=-PROV2
                      IEN(4,IE)=-PROV1
                      PROV1= - IEN(5,IE)
                      PROV2= - IEN(1,IE)
                      PROV3= - IEN(3,IE)
                      IEN(1,IE) = PROV2
                      IEN(3,IE) = PROV1
                      IEN(5,IE) = PROV3
                      LWRONG = .TRUE.
                      IEWRONG = IE
                      IEWRONGSUM = IEWRONGSUM + 1
                      WRITE(DBG%FHNDL,*) 'WRONG ELEMENT', IE
                      J1=INE(1,IE)
                      J2=INE(2,IE)
                      J3=INE(3,IE)
                      WRITE(DBG%FHNDL,*) 'NODENUMBERS I1=', J1
                      WRITE(DBG%FHNDL,*) 'NODENUMBERS I2=', J2
                      WRITE(DBG%FHNDL,*) 'NODENUMBERS I3=', J3
                      WRITE(DBG%FHNDL,*) 'XP1, YP1=', XP(J1), YP(J1)
                      WRITE(DBG%FHNDL,*) 'XP2, YP2=', XP(J2), YP(J2)
                      WRITE(DBG%FHNDL,*) 'XP3, YP3=', XP(J3), YP(J3)
                      WRITE(DBG%FHNDL,*) 'DXP1, DYP1=', DXP1, DYP1
                      WRITE(DBG%FHNDL,*) 'DXP2, DYP2=', DXP2, DYP2
                      WRITE(DBG%FHNDL,*) 'DXP3, DYP3=', DXP3, DYP3
                      WRITE(DBG%FHNDL,'(A40,6F15.8)') 'EDGELENGTHS OF THE WRONG ELEMENT', IEN(:,IE)
                   ELSE IF (TRIA(IE) .LT. THR) THEN 
                     write(DBG%FHNDL,*) 'IE=', IE, ' TRIA=', TRIA(IE)
                     write(DBG%FHNDL,*) 'DXP1=', DXP1, ' DXP3=', DXP3
                     write(DBG%FHNDL,*) 'DYP1=', DYP1, ' DYP3=', DYP3
                     write(DBG%FHNDL,*) 'I123=', I1, I2, I3
                     write(DBG%FHNDL,*) 'XP,YP(I1)=', XP(I1), YP(I1)
                     write(DBG%FHNDL,*) 'XP,YP(I2)=', XP(I2), YP(I2)
                     write(DBG%FHNDL,*) 'XP,YP(I3)=', XP(I3), YP(I3)
                     CALL WWM_ABORT('too small triangles')
                   END IF 

                   TL1 = SQRT(IEN(5,IE)**2 + IEN(6,IE)**2)
                   TL2 = SQRT(IEN(3,IE)**2 + IEN(4,IE)**2)
                   TL3 = SQRT(IEN(1,IE)**2 + IEN(2,IE)**2)
                   TMPTLMIN = MIN(TL1, TL2, TL3)
                   TMPTLMAX = MAX(TL1, TL2, TL3)

                   IF (TLMIN > TMPTLMIN) THEN
                      TLMIN = TMPTLMIN
                   END IF
                   IF (TLMAX < TMPTLMAX) THEN
                      TLMAX = TMPTLMAX
                   END IF

                   AVETA = AVETA+TRIA(IE)
                END DO

#ifdef MPI_PARALL_GRID
                CALL MPI_ALLREDUCE(TLMIN,TLMIN_GL,1,rtype,MPI_MIN,comm,ierr)
                CALL MPI_ALLREDUCE(TLMAX,TLMAX_GL,1,rtype,MPI_MAX,comm,ierr)
                CALL MPI_ALLREDUCE(AVETA,AVETA_GL,1,rtype,MPI_SUM,comm,ierr)
                AVETA = AVETA_GL / MyREAL(ne_global)
                AVETL = (TLMIN_GL+TLMAX_GL)/TWO
#else
                AVETA = AVETA/MyREAL(MNE)
                AVETL = (TLMIN+TLMAX)/TWO
#endif
#ifdef DEBUG
                WRITE(STAT%FHNDL,*) 'AVETL=', AVETL
                WRITE(STAT%FHNDL,*) 'TLMIN=', TLMIN
                WRITE(STAT%FHNDL,*) 'TLMAX=', TLMAX
                FLUSH(STAT%FHNDL)
#endif

                IF (LWRONG) THEN

#ifdef MPI_PARALL_GRID
                  IF (myrank == 0)  THEN
#endif
                    GRDCOR%FNAME = 'sysrenum.dat'
                    OPEN(GRDCOR%FHNDL, FILE=GRDCOR%FNAME, STATUS='UNKNOWN')
                    WRITE(GRDCOR%FHNDL,'(I10)') 0
                    WRITE(GRDCOR%FHNDL,'(I10)') MNP
                    IF (LSPHE) THEN
                      DO IP = 1, MNP
                        WRITE(GRDCOR%FHNDL,'(I10,3F15.8)') IP-1, XP(IP), YP(IP), DEP(IP)
                      END DO
                    ELSE
                      DO IP = 1, MNP
                        WRITE(GRDCOR%FHNDL,'(I10,3F15.6)') IP-1, XP(IP), YP(IP), DEP(IP)
                      END DO
                    ENDIF
                    WRITE(GRDCOR%FHNDL,'(I10)') MNE
                    DO IE = 1, MNE
                      WRITE(GRDCOR%FHNDL,'(5I10)') INE(1,IE)-1, INE(2,IE)-1, INE(3,IE)-1, 0, IE-1
                    END DO
                    CLOSE(GRDCOR%FHNDL)
#ifdef MPI_PARALL_GRID
                  END IF
#endif
                  WRITE(DBG%FHNDL,*) 'The Elements in your mesh are not correctly numbered!'
                  WRITE(DBG%FHNDL,*) 'New mesh is written to', TRIM(GRDCOR%FNAME)
                  WRITE(DBG%FHNDL,*) 'There are totally', IEWRONGSUM, 'wrong Elements' 
                  WRITE(DBG%FHNDL,*) 'The last wrong element has the number', IEWRONG
                  CALL WWM_ABORT('SPATIAL GRID - ELEMENTS NOT CORRECTLY NUMBERED')
                END IF
             CASE DEFAULT
               CALL WWM_ABORT('SPATIAL GRID - WRONG CASE - MUST BE 1D or 2D')
         END SELECT

         IF (LSPHE) THEN
            AVETL = AVETL*REARTH*PI/180.0
            DO IP = 1, MNP
               INVSPHTRANS(IP,1) = ONE/(DEGRAD*REARTH*COS(YP(IP)*DEGRAD))
               INVSPHTRANS(IP,2) = ONE/(DEGRAD*REARTH)
               IF ( ABS(COS(YP(IP))) .LT. THR ) THEN
                 CALL WWM_ABORT('SPATIAL GRID INVSPHETRANS IS SINGULAR')
               END IF
            END DO
         ELSE
           INVSPHTRANS = 1.
         END IF 

         IF (LTEST) THEN
            SELECT CASE (DIMMODE)
               CASE (1)
                  WRITE(STAT%FHNDL,101) AVETL, TLMIN, TLMAX
               CASE (2)
                  WRITE(STAT%FHNDL,102) AVETA, AVETL, TLMIN, TLMAX
                  IF (ITEST > 100) THEN
                     WRITE(STAT%FHNDL,*) ' The element area = ' 
                     DO IE = 1, MNE
                        WRITE(STAT%FHNDL,*) 'IE = ', IE, TRIA(IE)
                     END DO 
                  END IF
               CASE DEFAULT
            END SELECT
         END IF

         RETURN

101      FORMAT (1X,'The averege element length = ',F16.5/ &
     &           1X,'The minimum element length = ',F16.5/ &
     &           1X,'The maximum element length = ',F16.5/ )
102      FORMAT (1X,'The averege element area   = ',F16.5/ &
     &           1X,'The averege element length = ',F16.5/ &
     &           1X,'The minimum element length = ',F16.5/ &
     &           1X,'The maximum element length = ',F16.5/ )

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_TRIAD
      USE DATAPOOL
      IMPLICIT NONE
      TRI_ISP    = INT( LOG(TWO) / XISLN )
      TRI_ISP1   = TRI_ISP + 1
      TRI_WISP   = (2. - XIS**TRI_ISP) / (XIS**TRI_ISP1 - XIS**TRI_ISP)
      TRI_WISP1  = 1. - TRI_WISP
      TRI_ISM    = INT( LOG(0.5) / XISLN )
      TRI_ISM1   = TRI_ISM - 1
      TRI_WISM   = (XIS**TRI_ISM -0.5) / (XIS**TRI_ISM - XIS**TRI_ISM1)
      TRI_WISM1  = 1. - TRI_WISM
      TRI_ISBEGIN = MAX(1, 1-TRI_ISM1)
      IF (MESTR .eq. 1) THEN
        TRI_ARR(1)  = 0.1
        TRI_ARR(2)  = 2.2
        TRI_ARR(3)  = 10.
        TRI_ARR(4)  = 0.2
        TRI_ARR(5)  = 0.01
        IF (TRICO .GT. 0.)  TRI_ARR(1) = TRICO
        IF (TRIRA .GT. 0.)  TRI_ARR(2) = TRIRA
        IF (TRIURS .GT. 0.) TRI_ARR(5) = TRIURS
      END IF
      IF (MESTR .eq. 5) THEN
        TRI_ARR(1)  = 0.25
        TRI_ARR(2)  = 2.5
        TRI_ARR(3)  = 10.
        TRI_ARR(4)  = 0.2
        TRI_ARR(5)  = 0.01
        IF (TRICO .GT. 0.) TRI_ARR(1) = TRICO
        IF (TRIRA .GT. 0.) TRI_ARR(2) = TRIRA
        IF (TRIURS .GT. 0.) TRI_ARR(5) = TRIURS
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_SPECTRAL_GRID()
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER :: IS, ID
      INTEGER :: MSC1, MSC2
      REAL(rkind)    :: SSB, SPECTRAL_BANDWIDTH
      REAL(rkind)    :: TMP, CO1

      ALLOCATE( SPSIG(MSC), SPDIR(MDC), FR(MSC), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 5')
      SPSIG = zero
      SPDIR = zero
      FR    = zero

         ALLOCATE( COSTH(MDC), SINTH(MDC), COS2TH(MDC), SIN2TH(MDC), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 6')
         COSTH = zero
         SINTH = zero
         COS2TH = zero
         SIN2TH = zero

         ALLOCATE( SINCOSTH(MDC), SIGPOW(MSC,6), DS_BAND(0:MSC+1), DS_INCR(0:MSC+1), ID_NEXT(MDC), ID_PREV(MDC), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 7')
         SINCOSTH = zero
         SIGPOW = zero
         DS_BAND = zero
         DS_INCR = zero
         DO ID=1,MDC-1
           ID_NEXT(ID)=ID+1
         END DO
         ID_NEXT(MDC)=1
         DO ID=2,MDC
           ID_PREV(ID)=ID-1
         END DO
         ID_PREV(1)=MDC

         SGLOW  = PI2*FRLOW
         SGHIGH = PI2*FRHIGH

!2do check FRINTF for LOPTSIG

         IF (LOPTSIG) THEN
           SFAC   = 1.1_rkind
           FRINTF = 0.1_rkind
         ELSE 
           FRINTF = LOG(SGHIGH/SGLOW)/MyREAL(MSC-1) 
           SFAC   = EXP(FRINTF)
         END IF

         FRATIO = SFAC
         FRINTH = SQRT(SFAC)
         FR(1)  = FRLOW

         WRITE(STAT%FHNDL,*) 'RESOLUTION IN SIGMA SPACE AND FACTORS'
         WRITE(STAT%FHNDL,*) 'SGLOW', SGLOW
         WRITE(STAT%FHNDL,*) 'LOPTSIG', LOPTSIG
         WRITE(STAT%FHNDL,*) 'SFAC, FRINTF, FRINTH', SFAC, FRINTF, FRINTH

         DO IS = 2, MSC
           FR(IS) = FR(IS-1) * SFAC
         END DO

         SPSIG = FR * PI2 

         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'REL. FREQ. Distribution is =', FRINTF 

         IF ( ABS(FRINTF - .1)/FRINTF * 100. .GT. 1. ) THEN
           WRITE(DBG%FHNDL,*) 'Freq. resolution is not optimal for Snl4'
           WRITE(DBG%FHNDL,'(3F15.4)') 1. + FRINTF, ABS(FRINTF - .1)/FRINTF * 100.
           WRITE(DBG%FHNDL,*) 'rel. freq. res. should be 1.1 is now', 1. + FRINTF, 'ERROR IS:', ABS(FRINTF - .1)/FRINTF * 100.
         END IF  
         
         IF (MSC .GE. 2) THEN
           DS_BAND(0)     = SPSIG(2)- SPSIG(1)
           DS_BAND(1)     = DS_BAND(0)
           DS_BAND(MSC)   = SPSIG(MSC) - SPSIG(MSC-1)
           DS_BAND(MSC+1) = DS_BAND(MSC)
           DS_INCR(0)     = DS_BAND(0)
           DS_INCR(1)     = DS_BAND(0)
           DS_INCR(MSC)   = DS_BAND(MSC)
           DS_INCR(MSC+1) = DS_INCR(MSC)
           DO IS = 2, MSC-1 ! Bandwith at gridpoints
              DS_BAND(IS) = (SPSIG(IS)-SPSIG(IS-1))/2. + (SPSIG(IS+1)-SPSIG(IS))/2.
           END DO
           DO IS = 2, MSC ! Stepwidth between gridpoints K and K-1
              DS_INCR(IS) = SPSIG(IS) - SPSIG(IS-1)
           END DO
         END IF

         SPECTRAL_BANDWIDTH = (FRHIGH-FRLOW)
         SSB                =  (SPSIG(2)-SPSIG(1))/PI2
         MSCL               =  INT(SPECTRAL_BANDWIDTH/SSB)

         ALLOCATE(SPSIGL(MSCL));SPSIGL = ZERO
         SPSIGL(1) = FRLOW * PI2
         DO IS = 2, MSCL
           SPSIGL(IS) = SPSIGL(IS-1) + SSB
         ENDDO
!
!    *** the ratio of the consecutive frequency ... for quad
!
         IF ( MSC .GT. 3) THEN
           MSC2   = INT(MyREAL(MSC)/TWO)
           MSC1   = MSC2-1
           XIS    = SPSIG(MSC2)/SPSIG(MSC1)
           XISLN  = LOG(XIS)
           IF (SMETHOD .GT. 0 .AND. MESTR .GT. 0) CALL INIT_TRIAD
         ELSE
           IF (SMETHOD .GT. 0 .AND. MESNL .GT. 0) CALL WWM_ABORT('TOO LESS FREQ FOR SNL4 SET MESNL = 0')
           IF (SMETHOD .GT. 0 .AND. MESTR .GT. 0) CALL WWM_ABORT('TOO LESS FREQ FOR SNL3 SET MESTR = 0')
         END IF
!
!    *** frequency grid in [Hz]
!
         FR = SPSIG/PI2
! 
!    *** set the distribution of the dectional domain
! 
         IF (MDC == 1) THEN
           DDIR = 1._rkind
           SPDIR(1) = 0._rkind
         ELSE
            IF (LCIRD) THEN 
              DDIR = ABS(MAXDIR-MINDIR)/MyREAL(MDC)
              DO ID = 1, MDC
                 IF (LSTAG) THEN
                   SPDIR(ID) = MINDIR + DDIR * MyREAL(ID-1) + DDIR/2.0 
                 ELSE
                   SPDIR(ID) = MINDIR + DDIR * MyREAL(ID-1)
                 END IF
                 IF (SPDIR(ID) >= PI2) SPDIR(ID) = SPDIR(ID) - PI2
              END DO
            ELSE
              IF (LNAUTIN) THEN
                TMP = MAXDIR
                MAXDIR = MINDIR 
                MINDIR = TMP
              END IF 
              WRITE(STAT%FHNDL,*) 'MINDIR MAXDIR', MINDIR, MAXDIR, MINDIR*RADDEG, MAXDIR*RADDEG
              IF (MAXDIR.LT.MINDIR) MAXDIR = MAXDIR + PI2            
              DDIR = (MAXDIR-MINDIR) / MyREAL(MDC)                     
              DO ID = 1, MDC
                SPDIR(ID) = MINDIR + DDIR * MyREAL(ID-1)
              END DO
            END IF  
         END IF
!
!     *** set trig. in angular space 
!
         COSTH(:)    = COS(SPDIR(:))
         SINTH(:)    = SIN(SPDIR(:))
         COS2TH(:)   = COS(SPDIR(:))**2
         SIN2TH(:)   = SIN(SPDIR(:))**2
         SINCOSTH(:) = COS(SPDIR(:))*SIN(SPDIR(:))
!
!      *** set POWERS OF SPSIG
!
         SIGPOW(:,1) = SPSIG(:)
         SIGPOW(:,2) = SPSIG(:)**2
         SIGPOW(:,3) = SPSIG(:) * SIGPOW(:,2)
         SIGPOW(:,4) = SPSIG(:) * SIGPOW(:,3)
         SIGPOW(:,5) = SPSIG(:) * SIGPOW(:,4)
         SIGPOW(:,6) = SPSIG(:) * SIGPOW(:,5)
!
         FDIR = FRINTF * DDIR

         DELTH = PI2/MyREAL(MDC)

         ALLOCATE(DFIM(MSC), DFFR(MSC), DFFR2(MSC), DFIMOFR(MSC), FR5(MSC), FRM5(MSC), COFRM4(MSC), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_gridcf, allocate error 1')
         DFIM = ZERO
         DFFR = ZERO
         DFFR2 = ZERO 
         DFIMOFR = ZERO
         FR5 = ZERO
         FRM5 = ZERO 
         COFRM4 = ZERO

         if (LSOURCESWAM) then
           ALLOCATE(TH(MDC), stat=istat)
           th = zero
           DELTH = PI2/REAL(MDC)
           DELTR = DELTH*REARTH
           DO ID=1,MDC
             IF (LSTAG) THEN
               TH(ID) = REAL(ID-1)*DELTH
             ELSE
               TH(ID) = REAL(ID-1)*DELTH + 0.5*DELTH
             ENDIF
             COSTH(ID) = COS(TH(ID))
             SINTH(ID) = SIN(TH(ID))
           ENDDO
         endif

         CO1 = 0.5*(FRATIO-1.)*DELTH
         DFIM(1)= CO1*FR(1)
         DO IS=2,MSC-1
           DFIM(IS)=CO1 * (FR(IS)+FR(IS-1))
         ENDDO
         DFIM(IS)=CO1*FR(IS-1)

         DO IS = 1, MSC
           DFFR(IS)    = DFIM(IS)*FR(IS)
           DFFR2(IS)   = DFIM(IS)*FR(IS)**2
           DFIMOFR(IS) = DFIM(IS)/FR(IS)
           FR5(IS)     = FR(IS)**5
           FRM5(IS)    = ONE/FR5(IS)
           COFRM4(IS)  = COEF4*G9/FR(IS)**4
         END DO

         FLOGSPRDM1=1./LOG10(FRATIO)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEALLOC_SPECTRAL_GRID()
      USE DATAPOOL
      IMPLICIT NONE
      DEALLOCATE(DFIM, DFFR, DFFR2, DFIMOFR, FR5)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
