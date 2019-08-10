#include "wwm_functions.h"
#define DEBUG
#undef DEBUG
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SPECTRAL_SHAPE(SPPAR,ACLOC,LDEBUG,CALLFROM, OPTI)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(OUT)    ::  ACLOC(MSC,MDC)
      REAL(rkind), INTENT(INOUT)  ::  SPPAR(8)
      CHARACTER(LEN=*), INTENT(IN) :: CALLFROM
      LOGICAL, INTENT(IN) :: LDEBUG, OPTI
      IF (OPTI) THEN
        CALL OPTI_SPECTRAL_SHAPE(SPPAR,ACLOC,LDEBUG,CALLFROM)
      ELSE
        CALL KERNEL_SPECTRAL_SHAPE(SPPAR,ACLOC,LDEBUG,CALLFROM)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_ESTIMATE_PER_DIR_SHAPE(SPPAR, ACLOC, HS, TM, DM)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(IN)  :: SPPAR(8)
      REAL(rkind), INTENT(IN)  :: ACLOC(MSC,MDC)
      REAL(rkind), INTENT(OUT) :: HS, TM, DM
      REAL(rkind) :: DEPLOC, CURTXYLOC(2)
      REAL(rkind) :: WKLOC(MSC)
      REAL(rkind) :: FPP, TPP, CPP, WNPP, CGPP, KPP, LPP
      REAL(rkind) :: PEAKDSPR, PEAKDM, DPEAK, TPPD, KPPD, CGPD, CPPD
      REAL(rkind) :: TM01, TM02, TM10, KLM, WLM
      REAL(rkind) :: ETOTS, ETOTC, DSPR
      REAL(rkind) :: SPSIGLOC, WVN, WVC, WVK, WVCG
      integer ISMAX, IS
      DEPLOC=10000
      ISMAX=MSC
      CURTXYLOC=ZERO
      DO IS=1,MSC
        SPSIGLOC = SPSIG(IS)
        CALL WAVEKCG(DEPLOC,SPSIGLOC,WVN,WVC,WVK,WVCG)
        WKLOC(IS)=WVK
      END DO
      CALL MEAN_PARAMETER_LOC(ACLOC,CURTXYLOC,DEPLOC,WKLOC,ISMAX,HS,TM01,TM02,TM10,KLM,WLM)
      IF (SPPAR(5) .gt. 0) THEN
        CALL PEAK_PARAMETER_LOC(ACLOC,DEPLOC,ISMAX,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
        TM=TPP
        DM=PEAKDM
      ELSE
        CALL MEAN_DIRECTION_AND_SPREAD_LOC(ACLOC,ISMAX,ETOTS,ETOTC,DM,DSPR)
        TM=TM01
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OPTI_SPECTRAL_SHAPE(SPPAR,ACLOC,LDEBUG,CALLFROM)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(OUT)    ::  ACLOC(MSC,MDC)
      REAL(rkind), INTENT(INOUT)  ::  SPPAR(8)
      CHARACTER(LEN=*), INTENT(IN) :: CALLFROM
      LOGICAL, INTENT(IN) :: LDEBUG
      REAL(rkind) :: HS, TM, DM, TheErr, DeltaPer, Tper
      REAL(rkind) :: DiffAng, DEG, ADIR
      REAL(rkind) :: SPPARwork1(8), SPPARwork2(8), SPPARwork(8)
      integer :: iIter, nbIter, eSign, IS, ID
      REAL(rkind) :: eSum
      IF (ABS(SPPAR(5)) .eq. 3) THEN
        CALL KERNEL_SPECTRAL_SHAPE(SPPAR,ACLOC,LDEBUG,CALLFROM)
        RETURN
      END IF
      SPPARwork1=SPPAR
      SPPARwork2=SPPAR
      Tper=SPPAR(2)
      CALL KERNEL_SPECTRAL_SHAPE(SPPAR,ACLOC,LDEBUG,CALLFROM)
      CALL COMPUTE_ESTIMATE_PER_DIR_SHAPE(SPPAR, ACLOC, HS, TM, DM)
      DeltaPer=Tper - TM
      IF (TM < Tper) THEN
        eSign=1
      ELSE
        eSign=-1
      END IF
      SPPARwork=SPPAR
      SPPARwork(2)=SPPAR(2) + DeltaPer
      DO
        CALL KERNEL_SPECTRAL_SHAPE(SPPARwork,ACLOC,LDEBUG,CALLFROM)
        CALL COMPUTE_ESTIMATE_PER_DIR_SHAPE(SPPAR, ACLOC, HS, TM, DM)
        TheErr=(TM - Tper)*eSign
        IF (TheErr > 0) THEN
          EXIT
        END IF
        SPPARwork(2)=SPPARwork(2) + DeltaPer
      END DO
      IF (eSign .eq. 1) THEN
        SPPARwork1=SPPAR
        SPPARwork2=SPPARwork
      ELSE
        SPPARwork1=SPPARwork
        SPPARwork2=SPPAR
      END IF
      nbIter=20
      iIter=0
      DO
        SPPARwork=0.5_rkind*SPPARwork1 + 0.5_rkind*SPPARwork2
        CALL KERNEL_SPECTRAL_SHAPE(SPPARwork,ACLOC,LDEBUG,CALLFROM)
        CALL COMPUTE_ESTIMATE_PER_DIR_SHAPE(SPPAR, ACLOC, HS, TM, DM)
        IF (TM > Tper) THEN
          SPPARwork2=SPPARwork
        ELSE
          SPPARwork1=SPPARwork
        END IF
        iIter=iIter + 1
        IF (iIter > nbIter) THEN
          EXIT
        END IF
      END DO
      CALL KERNEL_SPECTRAL_SHAPE(SPPARwork,ACLOC,LDEBUG,CALLFROM)
      CALL COMPUTE_ESTIMATE_PER_DIR_SHAPE(SPPAR, ACLOC, HS, TM, DM)
      SPPARwork(1)=SPPAR(1)*(SPPAR(1)/HS)
      CALL KERNEL_SPECTRAL_SHAPE(SPPARwork,ACLOC,LDEBUG,CALLFROM)
      DO IS=1,MSC
        eSum=sum(ACLOC(IS,:))
      END DO
      CALL DEG2NAUT (SPPAR(3), DEG, LNAUTIN)
      ADIR = DEG * DEGRAD
      DO ID=1,MDC
        eSum=sum(ACLOC(:,ID))
        DiffAng=(360.0_rkind/PI2)*(SPDIR(ID) - ADIR)
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE KERNEL_SPECTRAL_SHAPE(SPPAR,ACLOC,LDEBUG,CALLFROM)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(OUT)    ::  ACLOC(MSC,MDC)
      REAL(rkind), INTENT(INOUT)  ::  SPPAR(8)
      CHARACTER(LEN=*), INTENT(IN) :: CALLFROM
      LOGICAL, INTENT(IN) :: LDEBUG

      INTEGER  ID, IS, LSHAPE, ITPER, ISP, ISIGMP
      REAL(rkind) ::   APSHAP, AUX1, AUX2, AUX3, AM0, AM1, AS2, AS3, ETOT, VEC2DEG
      REAL(rkind) ::  COEFF, SYF , MPER, PKPER, DIFPER, EHFR
      REAL(rkind) ::  MS, DEG, ETOTS, ETOTC, FF, CPSHAP, PPSHAP, DM, EAD, DS
      REAL(rkind) ::  RA, SALPHA, SF, SF4, SF5, FPK, FPK4, EFTAIL, CDIR
      REAL(rkind) ::  GAMMA_FUNC, DSPR, AACOS, ADIR, EPTAIL, APTAIL, PPTAIL
      REAL(rkind) ::  OMEG, EFTOT, ETOTT, CTOT, TM1, TPEAK
      LOGICAL  LOGPM, LINCOUT
!     SPPARM(1), WBHS: Hs, sign. wave height
!     SPPARM(2), WBTP: Wave period given by the user (either peak or mean)
!     SPPARM(3), WBDM: average direction
!     SPPARM(4), WBDS: directional spread
!     SPPARM(5), WBSS: spectral shape (1-4), (1 - Pierson-Moskowitz, 2 - JONSWAP, 3 - BIN, 4 - Gauss) peak (+) or mean frequency (-)
!     SPPARM(6), WBDSMS: directional spreading in degree (2) or exponent (1)
!     SPPARM(7), WBGAUSS: gaussian width for the gauss spectrum 0.1
!     SPPARM(8), WBPKEN: peak enhancement factor for the JONSWAP spectra 3.3

      IF (LDEBUG) THEN
        WRITE(DBG%FHNDL,*) 'HS    PER    DIR    DPSR    SHAPE   DEGEXP    GAUSS   PEAK'
        WRITE(DBG%FHNDL,'(8F10.4)') SPPAR(8)
      ENDIF

      ETOT = 0.
      EFTOT = 0.
      ACLOC = 0.

      IF (SPPAR(1) .LT. THR .OR. SPPAR(2) .LT. THR .OR. SPPAR(4) .LT. THR) THEN
        ACLOC = 0.
        RETURN
      END IF

      IF (SPPAR(5) .LT. 0) THEN
        LSHAPE = -INT(SPPAR(5))
        LOGPM  = .FALSE.
      ELSE
        LSHAPE =  INT(SPPAR(5))
        LOGPM  = .TRUE.
      ENDIF
!
      PKPER = SPPAR(2)
      ITPER = 0

      IF (LSHAPE.EQ.3) THEN
!       select bin closest to given period
        DIFPER = 1.E10
        DO IS = 1, MSC
          IF (ABS(PKPER - PI2/SPSIG(IS)) .LT. DIFPER) THEN
            ISP = IS
            DIFPER = ABS(PKPER - PI2/SPSIG(IS))
          END IF
        ENDDO
      ENDIF
!
100   FPK  = (1./PKPER)
      FPK4 = FPK**4

      IF (LSHAPE.EQ.1) THEN
        SALPHA = SPPAR(1)**2*FPK4 * 5./16.
      ELSE IF (LSHAPE.EQ.2) THEN
        SALPHA = (SPPAR(1)**2 * FPK4) / ((0.06533*(SPPAR(8)**0.8015)+0.13467)*16.)
      ELSE IF (LSHAPE.EQ.4) THEN
        AUX1 = SPPAR(1)**2 / ( 16.* SQRT (PI2) * SPPAR(7))
        AUX3 = 2._rkind * SPPAR(7)**2
      ENDIF
!
      DO IS = 1, MSC
!
        IF (LSHAPE.EQ.1) THEN

          SF = SPSIG(IS) / PI2
          SF4 = SF**4
          SF5 = SF**5
          RA = (SALPHA/SF5)*EXP(-(5.*FPK4)/(4.*SF4))/(PI2*SPSIG(IS))
          ACLOC(IS,MDC) = RA

        ELSE IF (LSHAPE.EQ.2) THEN

          SF = SPSIG(IS)/(PI2)
          SF4 = SF**4
          SF5 = SF**5
          CPSHAP = 1.25_rkind * FPK4 / SF4

          IF (CPSHAP.GT.10._rkind) THEN
            RA = 0._rkind
          ELSE
            RA = (SALPHA/SF5) * EXP(-CPSHAP)
          ENDIF

          IF (SF .LT. FPK) THEN
            COEFF = 0.07_rkind
          ELSE
            COEFF = 0.09_rkind
          ENDIF

          IF (COEFF*FPK .GT. SMALL) THEN
             APSHAP =  0.5_rkind * ((SF-FPK) / (COEFF*FPK))**2
          ELSE
             APSHAP =  ZERO
          ENDIF

          IF (APSHAP .GT. 10._rkind) THEN
            SYF = 1.
          ELSE
            PPSHAP = EXP(-APSHAP)
            SYF = SPPAR(8)**PPSHAP
          ENDIF

          RA = SYF*RA/(SPSIG(IS)*PI2)

          ACLOC(IS,MDC) = RA

          IF (LDEBUG) WRITE(DBG%FHNDL,*) 'IS LOOP', IS, SF, FPK, SYF, RA
!
        ELSE IF (LSHAPE .EQ. 3) THEN

          IF (IS.EQ.ISP) THEN
            ISBIN = ISP
            ACLOC(IS,MDC) = ( SPPAR(1)**2 ) / ( 16. * SIGPOW(IS,2) * FRINTF )
          ELSE
            ACLOC(IS,MDC) = 0.
          END IF

        ELSE IF (LSHAPE .EQ. 4) THEN

          AUX2 = ( SPSIG(IS) - ( PI2 / PKPER ) )**2
          RA = AUX1 * EXP ( -1. * AUX2 / AUX3 ) / SPSIG(IS)
          ACLOC(IS,MDC) = RA

        ELSE
          CALL WWM_ABORT('Wrong type for frequency shape 1 - 4')
        ENDIF

      END DO
      MPER = 0.
      IF (.NOT.LOGPM.AND.ITPER.LT.100) THEN
        ITPER = ITPER + 1
!       calculate average frequency
        AM0 = 0.
        AM1 = 0.
        DO IS = 1, MSC
          AS2 = ACLOC(IS,MDC) * SIGPOW(IS,2)
          AS3 = AS2 * SPSIG(IS)
          AM0 = AM0 + AS2
          AM1 = AM1 + AS3
        ENDDO
!       contribution of tail to total energy density
        PPTAIL = PTAIL(1) - 1.
        APTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
        AM0 = AM0 * FRINTF + APTAIL * AS2
        PPTAIL = PTAIL(1) - 2.
        EPTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
        AM1 = AM1 * FRINTF + EPTAIL * AS3
!       Mean period:
        IF ( AM1.GT.THR) THEN
          MPER = PI2 * AM0 / AM1
        ELSE
          MPER = ZERO 
        ENDIF
!        write(*,'(I10,4F15.8)') ITPER, MPER, &
!     &            ABS(MPER-SPPAR(2)) .GT. 0.01*SPPAR(2), &
!     &            (SPPAR(2) / MPER) * PKPER
        IF (ABS(MPER-SPPAR(2)) .GT. 0.01*SPPAR(2) .AND. MPER .GT. THR) THEN
!         modification suggested by Mauro Sclavo
          PKPER = (SPPAR(2)/MPER) * PKPER
          GOTO 100
        ELSE
          PKPER = ZERO 
        ENDIF 
      ELSE IF (ITPER.GE.100) THEN
        WRITE (STAT%FHNDL,*) 'No convergence calculating the spectrum'
        FLUSH(STAT%FHNDL)
      ENDIF

      CALL DEG2NAUT (SPPAR(3), DEG, LNAUTIN)

      ADIR = DEG * DEGRAD

      IF (INT(SPPAR(6)) .EQ. 1) THEN
        DSPR = PI * SPPAR(4) / 180._rkind
        MS = MAX (DSPR**(-2) - TWO, 1._rkind)
      ELSE
        MS = SPPAR(4)
      ENDIF

      IF (MS.LT.12._rkind) THEN
        CTOT = (2._rkind**MS)*(GAMMA_FUNC(0.5_rkind*MS+1._rkind))**2/(PI*GAMMA_FUNC(MS+1._rkind))
      ELSE
        CTOT =  SQRT(0.5_rkind*MS/PI)/(1._rkind-0.25_rkind/MS)
      ENDIF
      LINCOUT = .FALSE.
      DO ID = 1, MDC
        AACOS = COS(SPDIR(ID) - ADIR)
        !write(*,*) aacos, spdir(id), adir
        IF (AACOS .GT. 0._rkind) THEN
          CDIR = CTOT * MAX (AACOS**MS, THR)
          IF (.NOT.LCIRD) THEN
            IF (AACOS .GE. COS(DDIR)) THEN
              LINCOUT = .TRUE.
              !WRITE(*,*) 'ERROR', AACOS, COS(DDIR) 
              !STOP 'AVERAGE DIRECTION IS OUT OF SECTOR'
            ENDIF
          ENDIF
        ELSE
          CDIR = 0._rkind
        ENDIF
        DO IS = 1, MSC
          ACLOC(IS,ID) = CDIR * ACLOC(IS,MDC)
          !write(*,'(2I10,2F15.8)') is, id, cdir, ACLOC(IS,MDC)
        ENDDO
      ENDDO
!
! Finished creating
!
! Integral Parameters of the Input Spectra ...
! AR: 2DO: Optimize the output routine for this task ...
!
      IF (LDEBUG) THEN

        ETOT = 0.
        ETOTC = 0.
        ETOTS = 0.
 
        EFTAIL = 1.0 / (PTAIL(1)-1.0)

        IF (MSC .GE. 2) THEN
          DO ID = 1, MDC
            DO IS = 2, MSC
              DS = SPSIG(IS) - SPSIG(IS-1)
              EAD = 0.5*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS*DDIR
              ETOT = ETOT + EAD
              ETOTC  = ETOTC + EAD * COS(SPDIR(ID))
              ETOTS  = ETOTS + EAD * SIN(SPDIR(ID))
            END DO
            IF (MSC > 3) THEN
              EHFR = ACLOC(MSC,ID) * SPSIG(MSC)
              ETOT = ETOT + DDIR * EHFR * SPSIG(MSC) * EFTAIL
            ENDIF
          END DO
        ELSE
          DS = SGHIGH - SGLOW
          DO ID = 1, MDC
            EAD = ACLOC(1,ID) * DS * DDIR
            ETOT = ETOT + EAD
            ETOTC  = ETOTC + EAD * COS(SPDIR(ID))
            ETOTS  = ETOTS + EAD * SIN(SPDIR(ID))
          END DO
        END IF

        IF (ETOT > THR ) THEN
           DM    = VEC2DEG (ETOTC, ETOTS)
           CALL DEG2NAUT(DM,DEG,LNAUTOUT)
           DM = DEG
           FF = MIN (ONE, SQRT(ETOTC*ETOTC+ETOTS*ETOTS)/ETOT)
           DSPR = SQRT(2.-2.*FF) * 180./PI
        ELSE
          DM = 0.
          DSPR = 0.
        END IF

        ETOTT = 0.0
        EFTOT = 0.0
        EFTAIL = PTAIL(3)

        DO ID = 1, MDC
          DO IS = 1, MSC
            OMEG = SPSIG(IS)
            EAD = FRINTF * SIGPOW(IS,2) * ACLOC(IS,ID)
            ETOTT = ETOTT + EAD
            EFTOT = EFTOT + EAD * OMEG
          END DO
        END DO
        IF (EFTOT > VERYSMALL) THEN
          TM1 = PI2*ETOTT/EFTOT
        ELSE
          TM1 = 0.0
        END IF

        ETOTT = 0.0
        ISIGMP = -1
        DO IS = 1, MSC
         EAD = 0.0
         DO ID = 1, MDC
            EAD = EAD + SPSIG(IS)*ACLOC(IS,ID)*DDIR
         ENDDO
         IF (EAD > ETOTT) THEN
            ETOTT = EAD
            ISIGMP = IS
          END IF
        END DO
        IF (ISIGMP > 0) THEN
          TPEAK = 1./(SPSIG(ISIGMP)/PI2)
        ELSE
          TPEAK = 0.0
        END IF

        WRITE (STAT%FHNDL,'("+TRACE...",A)') 'GIVEN BOUNDARY SPECTRA AND RECALCULATED WAVE PARAMETERS'
        WRITE (STAT%FHNDL,'("+TRACE...",A)') 'THE DIFFERENCE IS MOSTLY DUE TO COARSE RESOLUTION IN SPECTRAL SPACE'
        WRITE (STAT%FHNDL,*) 'GIVEN     ', 'HS =', SPPAR(1)  
        WRITE (STAT%FHNDL,*) 'GIVEN     ', 'DM =', SPPAR(3)
        WRITE (STAT%FHNDL,*) 'GIVEN     ', 'DSPR =', SPPAR(4)
        WRITE (STAT%FHNDL,*) 'GIVEN     ', 'TM or TP', SPPAR(2)
        WRITE (STAT%FHNDL,*) 'SIMUL     ', 'HS =', 4. * SQRT(ETOT)
        WRITE (STAT%FHNDL,*) 'SIMUL     ', 'DM =',       DM 
        WRITE (STAT%FHNDL,*) 'SIMUL     ', 'DSPR =', DSPR
        WRITE (STAT%FHNDL,*) 'SIMUL     ', 'TM=', TM1, 'TPEAK=', TPEAK
        WRITE (STAT%FHNDL,*) 'TOT AC   =', SUM(ACLOC)
        WRITE (STAT%FHNDL,*) SPPAR
        FLUSH(STAT%FHNDL)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SPECTRUM_INT(ACLOC)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(OUT)   :: ACLOC(MSC,MDC)
      REAL(rkind)                :: MS(MSC), MS1, ADIR1, DS, EAD
      REAL(rkind)                :: INSPF(WBMSC)
      REAL(rkind)                :: INSPE(WBMSC)
      REAL(rkind)                :: INDIR(WBMSC)
      REAL(rkind)                :: INSPRD(WBMSC)
      REAL(rkind)                :: INMS(WBMSC)
      REAL(rkind)                :: SPCDIR(MSC)
      INTEGER                    :: IS, IS2, ID
      REAL(rkind)                :: CTOT(MSC), CDIRT, CDIR(MDC), CTOT1, CDIR1
      REAL(rkind)                :: DDACOS, DEG, DX, DIFFDX, YINTER
      REAL(rkind)                :: GAMMA_FUNC, ETOT, TM2
      REAL(rkind)                :: EFTOT, TM1, OMEG, PPTAIL, OMEG2
      REAL(rkind)                :: RA, ETAIL, EFTAIL
      REAL(rkind)                :: THD(MDC)
      REAL(rkind)                :: DTHD, RTH0 
      MS1 = WBDS
      ADIR1 = WBDM
      PSHAP(1) = 3.3
      PSHAP(2) = 0.07
      PSHAP(3) = 0.09
      PSHAP(4) = 0.01
      PSHAP(5) = 2.0
      PSHAP(6) = 0.01
!
!        2do Convert ... set the correct units ...
!
      DO IS = 1, WBMSC
        INSPF(IS)  = SFRQ(IS,1) * PI2
        INDIR(IS)  = SDIR(IS,1)
        INSPRD(IS) = SPRD(IS,1) * DEGRAD
        INSPE(IS)  = SPEG(IS,1,1) / PI2
      END DO
      ETOT = 0.0
      DO IS = 2, WBMSC
        DS = (INSPF(IS) - INSPF(IS-1)) 
        EAD = 0.5*(INSPE(IS) + INSPE(IS-1))*DS
        ETOT = ETOT + EAD
      END DO
      WRITE (STAT%FHNDL,*) 'HS - INPUTSPECTRA - 1', 4.0*SQRT(ETOT)
      ACLOC = 0.
      CALL INTERLIN (WBMSC, MSC, INSPF, SPSIG, INSPE, ACLOC(:,1))
      DO IS = 1, MSC
        IF (SPSIG(IS) .GT. INSPF(WBMSC)) THEN
           WRITE (STAT%FHNDL,*) 'Discrete Frequency is bigger then measured set FRMAX =', INSPF(WBMSC)/PI2
           WRITE (STAT%FHNDL,*) 'Setting all Action above the max. measured freq. zero'
           ACLOC(IS,1) = 0.0
        END IF
      END DO
      ETOT = 0.0
      DO IS = 2, MSC
        DS   = SPSIG(IS) - SPSIG(IS-1)
        EAD  = 0.5*(ACLOC(IS,1) + ACLOC(IS-1,1))*DS
        ETOT = ETOT + EAD
      END DO
      WRITE (STAT%FHNDL,*) 'HS - INPUTSPECTRA - INTERPOLATED', 4.0*SQRT(ETOT)
!
!        Convert to Wave Action if nessasary
!
      IF (LWBAC2EN) THEN
        DO IS = 1, MSC
          ACLOC(IS,1) = ACLOC(IS,1) / SPSIG(IS)
        END DO
      END IF
      ETOT = 0.0
      DO IS = 2, MSC
        DS   = SPSIG(IS) - SPSIG(IS-1)
        EAD  = 0.5 * (SPSIG(IS)*ACLOC(IS,1)+SPSIG(IS-1)*ACLOC(IS-1,1))*DS
        ETOT = ETOT + EAD
      END DO
      WRITE (STAT%FHNDL,*) 'HS - INPUTSPECTRA - WAVE ACTION', 4.0*SQRT(ETOT)
!
!        Convert from nautical to cartesian direction if nessasary and from deg2rad
!
      DO IS = 1, WBMSC
        CALL DEG2NAUT (INDIR(IS), DEG, LNAUTIN)
        INDIR(IS) = DEG
      END DO
!
!        Interpolate Directions in frequency space
!
      DO IS = 1, MSC
        DO IS2 = 1, WBMSC - 1
          IF (SPSIG(IS) .GT. INSPF(IS2) .AND. SPSIG(IS) .LT. INSPF(IS2+1)) THEN
            DX     = INSPF(IS2+1) - INSPF(IS2)
            DIFFDX = SPSIG(IS)    - INSPF(IS2)
            CALL INTERDIR( INDIR(IS2), INDIR(IS2+1), DX, DIFFDX, YINTER)
            SPCDIR(IS) = YINTER
            IF (SPSIG(IS) .GT. INSPF(WBMSC) ) SPCDIR(IS) = 0.0
          END IF
        END DO
        IF (SPSIG(IS) .GT. INSPF(WBMSC) ) SPCDIR(IS) = 0.0
      END DO
      CALL INTERLIN (WBMSC, MSC, INSPF, SPSIG, INDIR, SPCDIR)
      DO IS = 1, MSC
        IF ( SPSIG(IS) .GT. INSPF(WBMSC) ) SPCDIR(IS) = 0.
      END DO
      DO IS = 1, MSC
        DEG = SPCDIR(IS) * DEGRAD
        SPCDIR(IS) = DEG
      END DO
      IF (LINDSPRDEG) THEN
        DO IS = 1, WBMSC
          INMS(IS) = MAX (INSPRD(IS)**(-2) - TWO, ONE)
        END DO
      ELSE
        DO IS = 1, WBMSC
          INMS(IS) = INSPRD(IS)
        END DO
      END IF
!
!       Interpolate MS in Frequency Space, if LCUBIC than Cubic Spline Interpolation is used
!
      MS = 0.
      CALL INTERLIN (WBMSC, MSC, INSPF, SPSIG, INMS, MS)
      DO IS = 1, MSC
        IF ( SPSIG(IS) .GT. INSPF(WBMSC) ) MS(IS) = MS(IS-1)
      END DO
!
!        Construction of a JONSWAP Spectra
!
      IF (LPARMDIR) THEN
        IF (MS1.GT.10.) THEN
          CTOT1 = SQRT(MS1/(2.*PI)) * (1. + 0.25/MS1)
        ELSE
          CTOT1 = 2.**MS1 * (GAMMA_FUNC(1.+0.5*MS1))**2 / (PI * GAMMA_FUNC(1.+MS1))
        ENDIF
        DO IS = 1, MSC
          RA = ACLOC(IS,1)
          CDIRT = 0.
          DO ID = 1, MDC
            DDACOS = COS(SPDIR(ID) - ADIR1)
            IF (DDACOS .GT. 0.) THEN
              CDIR1 = CTOT1 * MAX (DDACOS**MS1, THR)
            ELSE
              CDIR1 = 0.
            ENDIF
            ACLOC(IS,ID) = RA * CDIR1
            CDIRT        = CDIRT + CDIR1 * DDIR
          END DO
        END DO
      ELSE
        DO IS = 1, MSC
          IF (MS(IS).GT.10.) THEN
            CTOT(IS) = SQRT(MS(IS)/(2.*PI)) * (1. + 0.25/MS(IS))
          ELSE
            CTOT(IS) = 2.**MS(IS) * (GAMMA_FUNC(1.+0.5*MS(IS)))**2.0 / (PI * GAMMA_FUNC(1.+MS(IS)))
          ENDIF
        END DO
        DO IS = 1, MSC
          RA = ACLOC(IS,1)
          CDIRT = 0.
          DO ID = 1, MDC
            DDACOS = COS(SPDIR(ID) - SPCDIR(IS))
            IF (DDACOS .GT. 0.) THEN
              CDIR(ID) = CTOT(IS) * MAX (DDACOS**MS(IS), ZERO)
            ELSE
              CDIR(ID) = 0.
            ENDIF
            ACLOC(IS,ID) = RA * CDIR(ID)
            CDIRT        = CDIRT + CDIR(ID) * DDIR
          ENDDO

          IF (100. - 1./CDIRT*100. .GT. 1.) THEN
            WRITE (STAT%FHNDL,*) 100 - 1./CDIRT*100., 'ERROR BIGGER THAN 1% IN THE DIRECTIONAL DISTTRIBUTION'
            WRITE (STAT%FHNDL,*) 'PLEASE CHECK THE AMOUNG OF DIRECTIONAL BINS AND THE PRESCRIBED DIRECTIONAL SPREADING'
          END IF
        ENDDO
      END IF
      ETOT = 0.
      DO ID = 1, MDC
        DO IS = 2, MSC
          DS = SPSIG(IS) - SPSIG(IS-1)
          EAD = 0.5*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS*DDIR
          ETOT = ETOT + EAD
        END DO
      END DO

      WRITE (STAT%FHNDL,*) 'HS - INPUTSPECTRA - AFTER 2D', 4.0*SQRT(ETOT)

      ETOT = 0.
      EFTOT = 0.
      PPTAIL = PTAIL(1) - 1.
      ETAIL  = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
      PPTAIL = PTAIL(1) - 2.
      EFTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))

      DO ID=1, MDC
         DO IS = 1, MSC
           OMEG = SPSIG(IS)
           EAD = FRINTF * SIGPOW(IS,2) * ACLOC(IS,ID)
           ETOT = ETOT + EAD
           EFTOT = EFTOT + EAD * OMEG
         ENDDO
         IF (MSC .GT. 3) THEN
           EAD = SIGPOW(MSC,2) * ACLOC(MSC,ID)
           ETOT = ETOT + ETAIL * EAD
           EFTOT = EFTOT + EFTAIL * OMEG * EAD
!           WRITE (*,*)  ETAIL * EAD, EFTAIL * OMEG * EAD, ACLOC(MSC,ID)
         ENDIF
      ENDDO
      IF (EFTOT.GT.THR) THEN
         TM1 = PI2 * ETOT / EFTOT
      ELSE
         TM1 = 0.
      ENDIF

      ETOT  = 0.
      EFTOT = 0.
      PPTAIL = PTAIL(1) - 1.
      ETAIL  = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
      PPTAIL = PTAIL(1) - 3.
      EFTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
      DO ID=1, MDC
         DO IS=1,MSC
           EAD  = SIGPOW(IS,2) * ACLOC(IS,ID) * FRINTF
           OMEG2 = SIGPOW(IS,2)
           ETOT  = ETOT + EAD
           EFTOT = EFTOT + EAD * OMEG2
         ENDDO
         IF (MSC .GT. 3) THEN
           EAD  = SIGPOW(MSC,2) * ACLOC(MSC,ID)
           ETOT  = ETOT  + ETAIL * EAD
           EFTOT = EFTOT + EFTAIL * OMEG2 * EAD
         ENDIF
      ENDDO
      IF (EFTOT .GT. THR) THEN
         TM2 = PI2 * SQRT(ETOT/EFTOT)
      ELSE
         TM2 = 0.
      END IF

      ETOT = 0.
      DO ID = 1, MDC
        DO IS = 2, MSC
          DS = SPSIG(IS) - SPSIG(IS-1)
          EAD = 0.5*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS*DDIR
          ETOT = ETOT + EAD
        END DO
      END DO

      WRITE (STAT%FHNDL,*) 'HS - INPUTSPECTRA - AFTER 2D', 4.0*SQRT(ETOT)
      WRITE (STAT%FHNDL,*) 'TM01, TM02 & HS', TM1, TM2, 4.0*SQRT(ETOT)

      !FLUSH(DBG%FHNDL)
      !FLUSH(STAT%FHNDL)

      IF (.FALSE.) THEN ! Write WW3 spectra of the input boundary condition ...

        DTHD=360./MDC
        RTH0=SPDIR(1)/DDIR
        DO ID = 1, MDC
          THD(ID)=DTHD*(RTH0+MyREAL(ID-1))
        END DO
        WRITE (4001,1944) 'WAVEWATCH III SPECTRA', MSC, MDC, 1, 'LAI ET AL'
        WRITE (4001,1945) (SPSIG(IS)*INVPI2,IS=1,MSC)
        WRITE (4001,1946) (MOD(2.5*PI-SPDIR(ID),PI2),ID=1,MDC)
        WRITE (4001,901) 'LAI SPEC', 0., 0., 0., 0., 0., 0., 0.
        WRITE (4001,902) ((ACLOC(IS,ID)*SPSIG(IS)/PI2*RHOW*G9,IS=1,MSC),ID=1,MDC)

      END IF

  901 FORMAT ('''',A10,'''',2F7.2,F10.1,2(F7.2,F6.1))
  902 FORMAT (7E11.3)
 1943 FORMAT ( '      File name : ',A,' (',A,')')
 1944 FORMAT ('''',A,'''',1X,3I6,1X,'''',A,'''')
 1945 FORMAT (8E10.3)
 1946 FORMAT (7E11.3)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_WAVE_BOUNDARY
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER :: IP, eIdx
      IF (LBCWA .OR. LBCSP) THEN ! Spectrum or parametric boundary condition
        IF (LINHOM) THEN
          DO IP = 1, IWBMNP
            eIdx = IWBNDLC(IP)
            AC2(:,:,eIdx) = WBAC(:,:,IP)
          END DO
        ELSE
          DO IP = 1, IWBMNP
            eIdx = IWBNDLC(IP)
            AC2(:,:,eIdx) = WBAC(:,:,1)
          END DO
        END IF
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_IFILE_IT(IFILE, IT)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: IFILE, IT
      REAL(rkind) :: DTMP
      INTEGER ITMP
      DTMP = (MAIN%TMJD-BND_TIME_ALL_FILES(1,1)) * DAY2SEC
      ITMP  = 0
      DO IFILE = 1, NUM_NETCDF_FILES_BND
        ITMP = ITMP + NDT_BND_FILE(IFILE)
        IF (ITMP .GT. INT(DTMP/SEBO%DELT)) EXIT
      END DO
      ITMP = SUM(NDT_BND_FILE(1:IFILE-1))
      IT   = NINT(DTMP/SEBO%DELT) - ITMP + 1
      IF (IT .GT. NDT_BND_FILE(IFILE)) THEN
        IFILE = IFILE + 1
        IT    = 1
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_BOUNDARY_CONDITION(WBACOUT)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(OUT)   :: WBACOUT(MSC,MDC,IWBMNP)
      INTEGER                    :: IP
!AR: WAVE BOUNDARY

!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by user (either peak or mean)
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!     SPPARM(5): spectral shape (1-4),
!         1 - Pierson-Moskowitz
!         2 - JONSWAP
!         3 - BIN
!         4 - Gauss
!         positive peak (+) or mean frequency (-)

!     SPPARM(6): directional spreading in degree (1) or exponent (2)
!     SPPARM(7): gaussian width for the gauss spectrum 0.1
!     SPPARM(8): peak enhancement factor for the JONSWAP spectra 3.
!
! Count number of active boundary points ...
!
      IF(LWW3GLOBALOUT) THEN
        IF (.NOT. ALLOCATED(WW3GLOBAL)) THEN
          ALLOCATE(WW3GLOBAL(8,MNP), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 7')
        END IF
      END IF
      IF (LBCWA) THEN ! Parametric Wave Boundary is prescribed
        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'Parametric Wave Boundary Condition is prescribed'
        IF (LINHOM) THEN ! Inhomogenous in space
          IF (LBCSE) THEN ! Unsteady in time
            SPPARM = 0.
            IF (IBOUNDFORMAT == 1) THEN  ! WWM
              CALL READWAVEPARWWM
            ELSE IF (IBOUNDFORMAT == 2) THEN ! FVCOM ... THIS WILL BE REPLACED BY SWAN TYPE BOUNDARY!
              CALL READWAVEPARFVCOM
            ELSE IF (IBOUNDFORMAT == 3) THEN ! WW3
#ifdef NCDF
              CALL READ_NETCDF_WW3_PARAM
#else
              CALL WWM_ABORT('compile with netcdf for IBOUNDFORMAT=3')
#endif
              CALL INTER_STRUCT_BOUNDARY(NDX_BND,NDY_BND,DX_BND,DY_BND,OFFSET_X_BND,OFFSET_Y_BND,SPPARM)
              IF (LWW3GLOBALOUT) CALL INTER_STRUCT_DOMAIN(NDX_BND,NDY_BND,DX_BND,DY_BND,OFFSET_X_BND,OFFSET_Y_BND,WW3GLOBAL)
            ELSE IF (IBOUNDFORMAT == 4) THEN ! WWM SPPARM netcdf file
!              Print *, 'Before READ_NETCDF_BOUNDARY_SPPARM'
#ifdef NCDF
              CALL READ_NETCDF_BOUNDARY_SPPARM
#else
              CALL WWM_ABORT('compile with netcdf for IBOUNDFORMAT=4')
#endif
            ELSE IF (IBOUNDFORMAT == 5) THEN ! WAM format of waves
               CALL WWM_ABORT('No possibility of using parametric boundary for IBOUNDFORMAT=5')
            END IF
          ELSE  ! Steady ...
            SPPARM = 0.
            IF (IBOUNDFORMAT == 1) THEN
              CALL READWAVEPARWWM
            ELSE IF (IBOUNDFORMAT == 2) THEN
              CALL READWAVEPARFVCOM
            ELSE IF (IBOUNDFORMAT == 3) THEN
#ifdef NCDF
              CALL READ_NETCDF_WW3_PARAM
#else
              CALL WWM_ABORT('compile with netcdf for IBOUNDFORMAT=3')
#endif
              CALL INTER_STRUCT_BOUNDARY(NDX_BND,NDY_BND,DX_BND,DY_BND,OFFSET_X_BND,OFFSET_Y_BND,SPPARM)
            END IF
          END IF ! LBCSE ...
          DO IP = 1, IWBMNP
            CALL SPECTRAL_SHAPE(SPPARM(:,IP),WBACOUT(:,:,IP),.FALSE.,'CALL FROM WB 1', .FALSE.)
          END DO
        ELSE ! Homogenous in space
          IF (IWBMNP .gt. 0) THEN
            IF (LBCSE) THEN ! Unsteady in time
              IF (IBOUNDFORMAT == 1) THEN
                CALL READWAVEPARWWM
              ELSE IF (IBOUNDFORMAT == 2) THEN
                CALL READWAVEPARFVCOM
              END IF
              CALL SPECTRAL_SHAPE(SPPARM(:,1),WBACOUT(:,:,1), .FALSE.,'CALL FROM WB 3', .FALSE.)
            ELSE ! Steady in time ...
              SPPARM = 0.
              IF (LMONO_IN) THEN
                SPPARM(1,1) = WBHS * SQRT(2.)
              ELSE
                SPPARM(1,1) = WBHS
              END IF
              SPPARM(2,1) = WBTP
              SPPARM(3,1) = WBDM
              SPPARM(4,1) = WBDS
              SPPARM(5,1) = WBSS
              SPPARM(6,1) = WBDSMS
              SPPARM(7,1) = WBGAUSS
              SPPARM(8,1) = WBPKEN
              CALL SPECTRAL_SHAPE(SPPARM(:,1),WBACOUT(:,:,1),.FALSE.,'CALL FROM WB 4', .TRUE.)
            END IF ! LBCSE
          END IF
        END IF ! LINHOM
      END IF
      IF (LBCSP) THEN ! Spectrum is prescribed
        IF (LINHOM) THEN ! The boundary conditions is not homogenous!
          IF (LBSP1D) CALL WWM_ABORT('No inhomogenous 1d spectra boundary cond. available') 
          IF (IBOUNDFORMAT == 1) THEN ! WWM
            !CALL READSPEC2D
            CALL WWM_ABORT('No inhomogenous 2d spectra boundary cond. available in WWM Format')
          END IF
          IF (IBOUNDFORMAT == 3) THEN ! WW3
            CALL GET_BINARY_WW3_SPECTRA(WBACOUT)
          END IF
          IF (IBOUNDFORMAT == 4) THEN ! WWM WBAC netcdf
#ifdef NCDF
            CALL READ_NETCDF_BOUNDARY_WBAC(WBACOUT)
#else
            CALL WWM_ABORT('compile with netcdf for IBOUNDFORMAT=4')
#endif
          END IF
          IF (IBOUNDFORMAT == 5) THEN ! WWM WBAC netcdf
            Print *, 'Before call to READ_GRIB_WAM_BOUNDARY_WBAC'
#ifdef GRIB_API_ECMWF
            CALL READ_GRIB_WAM_BOUNDARY_WBAC(WBACOUT)
#else
            CALL WWM_ABORT('compile with GRIB_API for IBOUNDFORMAT=5')
#endif
          END IF
          IF (IBOUNDFORMAT == 6) THEN ! WW3 netcdf
#ifdef NCDF
            CALL GET_NC_WW3_SPECTRA(WBACOUT)
#else
            CALL WWM_ABORT('compile with netcdf for IBOUNDFORMAT=4')
#endif
          END IF
        ELSE ! The boundary conditions is homogeneous in space !
          IF (LBSP1D) THEN ! 1-D Spectra is prescribed
            WRITE(STAT%FHNDL,'("+TRACE...",A)') '1d Spectra is given as Wave Boundary Condition'
            CALL READSPEC1D(LFIRSTREAD)
            CALL SPECTRUM_INT(WBACOUT(:,:,1))
          ELSE IF (LBSP2D) THEN ! 2-D Spectra is prescribed
            WRITE(STAT%FHNDL,'("+TRACE...",A)') '2d Spectra is given as Wave Boundary Condition'
            IF (IBOUNDFORMAT == 1) THEN
              CALL READSPEC2D
            ELSE IF (IBOUNDFORMAT == 3) THEN
              CALL GET_BINARY_WW3_SPECTRA(WBACOUT) 
            ELSE IF (IBOUNDFORMAT == 6) THEN
              CALL GET_NC_WW3_SPECTRA(WBACOUT) 
            END IF
            CALL SPECTRUM_INT(WBACOUT(:,:,1))
          END IF ! LBSP1D .OR. LBSP2D
        END IF ! LINHOM
      ENDIF ! LBCSP
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_WAVE_BOUNDARY_CONDITION
      USE DATAPOOL
      IMPLICIT NONE
      CHARACTER(len=29) :: CHR
      IF (LNANINFCHK) THEN
        WRITE(DBG%FHNDL,*) ' ENTERING SET_WAVE_BOUNDARY_CONDITION ',  SUM(AC2)
        IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN BOUNDARY CONDTITION l.1959')
      ENDIF
      IF ( MAIN%TMJD > SEBO%TMJD-1.E-8 .AND. MAIN%TMJD < SEBO%EMJD ) THEN ! Read next time step from boundary file ...
        IF (LBINTER) THEN
          CALL WAVE_BOUNDARY_CONDITION(WBACNEW)
          DSPEC   = (WBACNEW-WBACOLD)/SEBO%DELT*MAIN%DELT
          WBAC    =  WBACOLD
          WBACOLD =  WBACNEW
        ELSE ! .NOT. LBINTER
          CALL WAVE_BOUNDARY_CONDITION(WBAC)
        END IF ! LBINTER
        SEBO%TMJD = SEBO%TMJD + SEBO%DELT*SEC2DAY
      ELSE ! Interpolate in time ... no need to read ...
        IF (LBINTER) THEN
          WBAC = WBAC + DSPEC
        END IF
      END IF
      CALL SET_WAVE_BOUNDARY
      IF (LNANINFCHK) THEN
        WRITE(DBG%FHNDL,*) ' FINISHED WITH BOUNDARY CONDITION ',  SUM(AC2)
        IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN BOUNDARY CONDTITION l.1959')
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_WAVE_BOUNDARY_CONDITION
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind)            :: DTMP
      INTEGER                :: ITMP
      LOGICAL                :: DoAllocate
      WRITE(STAT%FHNDL, *) 'Begin of INIT_WAVE_BOUNDARY_CONDITION'
!TODO: Makes sure initial condition work also when no wave boundary is set ...
      IF (IBOUNDFORMAT == 3) THEN
        IF (LBCSP) THEN ! Spectrum is prescribed
          CALL INIT_BINARY_WW3_SPECTRA
        END IF
        IF (LBCWA) THEN ! Parametric Wave Boundary is prescribed
#ifdef NCDF
          CALL INIT_NETCDF_WW3_WAVEPARAMETER
#else
          CALL WWM_ABORT('Compile with netcdf For WW3 bdcons')
#endif
        END IF
      END IF
      IF (IBOUNDFORMAT == 4) THEN
#ifdef NCDF
        CALL INIT_NETCDF_BOUNDARY_WWM
#else
        CALL WWM_ABORT('Compile with netcdf for IBOUNDFORMAT=4')
#endif
      END IF
      IF (IBOUNDFORMAT == 5) THEN
#ifdef GRIB_API_ECMWF
        CALL INIT_GRIB_WAM_BOUNDARY
#else
        CALL WWM_ABORT('Compile with GRIB_API for IBOUNDFORMAT=5')
#endif
      END IF
      IF (IBOUNDFORMAT == 6) THEN
#ifdef NCDF
      IF (LBCSP) THEN ! Spectrum is prescribed
        CALL INIT_NC_WW3_SPECTRA
      ELSE
        CALL WWM_ABORT('IBOUNDFORMAT == 6 is for prescribed spectra only')
      END IF
#else
        CALL WWM_ABORT('Compile with netcdf for IBOUNDFORMAT=6')
#endif
      END IF
      IF ((IBOUNDFORMAT .eq. 4).or.(IBOUNDFORMAT .eq. 5).or.LEXPORT_BOUC_WW3.or.BOUC_NETCDF_OUT_PARAM.or.BOUC_NETCDF_OUT_SPECTRA) THEN
#ifdef MPI_PARALL_GRID
        CALL SETUP_BOUNDARY_SCATTER_REDUCE_ARRAY
#endif
#ifdef MPI_PARALL_GRID
        IF (myrank .eq. rank_boundary) THEN
          DoAllocate=.TRUE.
        ELSE
          DoAllocate=.FALSE.
        END IF
#else
        DoAllocate=.TRUE.
#endif
        WRITE(STAT%FHNDL, *) 'DoAllocate=', DoAllocate
        IF (DoAllocate) THEN
          IF (LBCWA .or. BOUC_NETCDF_OUT_PARAM) THEN
            IF (.NOT. ALLOCATED(SPPARM_GL)) THEN
              allocate(SPPARM_GL(8,IWBMNPGL), stat=istat)
              IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 26')
            END IF
          END IF 
          IF (LBCSP .or. LEXPORT_BOUC_WW3 .or. BOUC_NETCDF_OUT_SPECTRA) THEN
            IF (.NOT. ALLOCATED(WBAC_GL)) THEN
              allocate(WBAC_GL(MSC,MDC,IWBMNPGL), stat=istat)
              IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 27')
            END IF
          END IF
        END IF
      END IF
      WRITE(STAT%FHNDL, *) 'LEXPORT_BOUC_WW3=', LEXPORT_BOUC_WW3
      WRITE(STAT%FHNDL, *) 'IWBMNP=', IWBMNP
      WRITE(STAT%FHNDL, *) 'allocated(WBAC)=', allocated(WBAC)
      WRITE(STAT%FHNDL, *) 'allocated(WBAC_GL)=', allocated(WBAC_GL)
      IF (LBCWA .or. LBCSP) THEN
        CALL WAVE_BOUNDARY_CONDITION(WBAC)
        IF (LBINTER) WBACOLD = WBAC
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef NCDF
      SUBROUTINE INIT_NETCDF_WW3_WAVEPARAMETER
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE

      INTEGER :: IT, IFILE, IVAR, BND_NCID
      INTEGER :: ILON_ID, ILAT_ID, ITIME_ID, I, J, COUNTER
      REAL(rkind), ALLOCATABLE :: BND_TIME(:)
      character (len = *), parameter :: CallFct = "INIT_NETCDF_WW3_WAVEPARAMETER"
      integer, dimension(nf90_max_var_dims) :: dimIDs
      INTEGER MAX_NDT_BND_FILE, iProc, eSize
      INTEGER iVect(4)
      REAL(rkind) rVect(4)
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_BOUND .or. (myrank .eq. 0)) THEN
# endif
        CALL TEST_FILE_EXIST_DIE("Missing WW3 boundary file : ", TRIM(WAV%FNAME))
        OPEN(WAV%FHNDL,FILE=WAV%FNAME,STATUS='OLD')
!        WRITE(STAT%FHNDL,*) WAV%FHNDL, WAV%FNAME, BND%FHNDL, BND%FNAME

        NUM_NETCDF_FILES_BND = 0
        DO
          READ( WAV%FHNDL, *, IOSTAT = ISTAT )
          IF ( ISTAT /= 0 ) EXIT
          NUM_NETCDF_FILES_BND = NUM_NETCDF_FILES_BND + 1
        END DO
        REWIND(WAV%FHNDL)
!        WRITE(STAT%FHNDL,*) 'NUM_NETCDF_FILES_BND=', NUM_NETCDF_FILES_BND

        NUM_NETCDF_FILES_BND = NUM_NETCDF_FILES_BND / NUM_NETCDF_VAR_TYPES
        WRITE(STAT%FHNDL,*) 'NUM_NETCDF_FILES_BND=', NUM_NETCDF_FILES_BND
        WRITE(STAT%FHNDL,*) 'NUM_NETCDF_VAR_TYPES=', NUM_NETCDF_VAR_TYPES
        ALLOCATE(NETCDF_FILE_NAMES_BND(NUM_NETCDF_FILES_BND,NUM_NETCDF_VAR_TYPES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 9')
        DO IT = 1, NUM_NETCDF_FILES_BND
          DO IVAR = 1, NUM_NETCDF_VAR_TYPES
            READ( WAV%FHNDL, *) NETCDF_FILE_NAMES_BND(IT,IVAR)
          END DO
        END DO
        CLOSE (WAV%FHNDL)
!
! four files are read to set up the wave spectra Hs, Tm01, Dir, Sprd
!
        ALLOCATE(NDT_BND_FILE(NUM_NETCDF_FILES_BND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 10')
        NDT_BND_FILE = 0
!        DO IFILE = 1, NUM_NETCDF_FILES_BND
!          WRITE(STAT%FHNDL,'(I10,10X,5A30)') IFILE, NETCDF_FILE_NAMES_BND(IFILE,:)
!        END DO
!
! check number of time steps in netcdf file ... it is assumed that all files have the same ammount of time steps ...
!
        DO IFILE = 1, NUM_NETCDF_FILES_BND
          WRITE(STAT%FHNDL,*) ifile, TRIM(NETCDF_FILE_NAMES_BND(IFILE,1))
          FLUSH(STAT%FHNDL)
          CALL TEST_FILE_EXIST_DIE("Missing ww3 boundary condition file : ", TRIM(NETCDF_FILE_NAMES_BND(IFILE,1)))
          ISTAT = NF90_OPEN(TRIM(NETCDF_FILE_NAMES_BND(IFILE,1)), NF90_NOWRITE, BND_NCID)
          IF (ISTAT /= 0) THEN
            Print *, 'Error while trying to open netcdf file'
            Print *, 'FILE=', TRIM(NETCDF_FILE_NAMES_BND(IFILE,1))
            Print *, 'One possible error is that the file is NC4 but you linked to NC3'
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)
          END IF

          ISTAT = nf90_inq_varid(BND_NCID, 'time', ITIME_ID)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)

          ISTAT = NF90_INQUIRE_VARIABLE(BND_NCID, ITIME_ID, dimids = dimids)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, ISTAT)

          ISTAT = nf90_inquire_dimension(BND_NCID, dimIDs(1), len = NDT_BND_FILE(IFILE))
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, ISTAT)

!          WRITE(STAT%FHNDL,*) IFILE, NDT_BND_FILE(IFILE)
        END DO
!
! check dimensions in the netcdf ... again it is assumed that this is not changing for all files ...
!
        ISTAT = nf90_inq_varid(BND_NCID, 'longitude', ILON_ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(BND_NCID, ILON_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, ISTAT)

        ISTAT = nf90_inquire_dimension(BND_NCID, dimIDs(1), len = NDX_BND)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, ISTAT)

        ISTAT = nf90_inq_varid(BND_NCID, 'latitude', ILAT_ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(BND_NCID, ILAT_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, ISTAT)

        ISTAT = nf90_inquire_dimension(BND_NCID, dimIDs(1), len = NDY_BND)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, ISTAT)

!        WRITE(STAT%FHNDL,*) 'Number of Gridpoints', NDX_BND, NDY_BND

        ALLOCATE (COORD_BND_X(NDX_BND), COORD_BND_Y(NDY_BND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 11')
!
! read coordinates from files ....
!
        ISTAT = NF90_GET_VAR(BND_NCID, ILON_ID, COORD_BND_X)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)

        ISTAT = NF90_GET_VAR(BND_NCID, ILAT_ID, COORD_BND_Y)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, ISTAT)
!
! estimate offset ...
!
        OFFSET_X_BND = MINVAL(COORD_BND_X)
        OFFSET_Y_BND = MINVAL(COORD_BND_Y)
!
! resolution ...
!
        DX_BND  = ABS(MAXVAL(COORD_BND_X) - MINVAL(COORD_BND_X))/(NDX_BND-1)
        DY_BND  = ABS(MAXVAL(COORD_BND_Y) - MINVAL(COORD_BND_Y))/(NDY_BND-1)
!
! close netcdf file ...
!
        ISTAT = NF90_CLOSE(BND_NCID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, ISTAT)
!
! total number of time steps ... in all files
!
        NDT_BND_ALL_FILES = 0
!        write(STAT%FHNDL,*) NUM_NETCDF_FILES_BND
        DO IT = 1, NUM_NETCDF_FILES_BND
          NDT_BND_ALL_FILES = NDT_BND_ALL_FILES + NDT_BND_FILE(IT)
!          write(STAT%FHNDL,*) it, NDT_BND_FILE(it)
        END DO
!        WRITE(STAT%FHNDL,*) NDT_BND_ALL_FILES, NDT_BND_FILE

        WRITE(STAT%FHNDL, *) 'INIT_NETCDF_WW3_WAVEPARAMETER'
        MAX_NDT_BND_FILE=MAXVAL(NDT_BND_FILE)
        ALLOCATE (BND_TIME_ALL_FILES(NUM_NETCDF_FILES_BND,MAX_NDT_BND_FILE), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 12')
        BND_TIME_ALL_FILES = ZERO
!
! read all time steps in the proper format and transform in wwm time line
!
        DO IFILE = 1, NUM_NETCDF_FILES_BND
          CALL TEST_FILE_EXIST_DIE("Missing ww3 boundary condition file : ", TRIM(NETCDF_FILE_NAMES_BND(IFILE,1)))
          ISTAT = NF90_OPEN(NETCDF_FILE_NAMES_BND(IFILE,1),NF90_NOWRITE,BND_NCID)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, ISTAT)

          ALLOCATE (BND_TIME(NDT_BND_FILE(IFILE)), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 13')
          BND_TIME = ZERO
! MDS: It looks dangerous to use previous id.
          ISTAT = NF90_GET_VAR(BND_NCID,ITIME_ID,BND_TIME)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 15, ISTAT)

          DO IT = 1, NDT_BND_FILE(IFILE)
            BND_TIME_ALL_FILES(IFILE,IT) = BND_TIME(IT)
!             CALL CT2MJD('19000101.000000',DTMP1)
!             CALL CT2MJD('19900101.000000',DTMP2)
!             CALL MJD2CT(DTMP1,chrdate)
!             WRITE(*,*) '19000101.000000', DTMP1, chrdate
!             CALL MJD2CT(DTMP2,chrdate)
!             WRITE(*,*) '19900101.000000', DTMP2, chrdate
!             CALL MJD2CT(0.0_rkind,chrdate)
!             WRITE(*,*) '00000000.000000', 0.0_rkind, chrdate
!             WRITE(*,*) BND_TIME_ALL_FILES(1,1), DT_DIFF_19901900
!             IF (IT == 1 .AND. IFILE ==1) WRITE(*,*) DTMP1, DTMP2, DTMP1+DT_DIFF_19901900
!             IF (IT == 1 .AND. IFILE ==1) WRITE(*,*) IFILE, IT, BND_TIME(IT), chrdate
          END DO
          DEALLOCATE(BND_TIME)
        END DO

        BND_TIME_ALL_FILES = BND_TIME_ALL_FILES + DT_DIFF_19901900

        IF (LWRITE_ALL_WW3_RESULTS) THEN
          OPEN(3010, FILE  = 'sysglobalboundary.dat', STATUS = 'UNKNOWN')
          WRITE (3010, '(I10)') 0
          WRITE (3010, '(I10)') NDX_BND * NDY_BND
          COUNTER = 0
          DO I = 1, NDY_BND
            DO J = 1, NDX_BND
              WRITE (3010, '(I10,3F15.4)') COUNTER, OFFSET_X_BND+(J-1)*DX_BND,OFFSET_Y_BND+(I-1)*DY_BND, 0.0
              COUNTER = COUNTER + 1
            END DO
          END DO
          WRITE (3010, *) (NDX_BND-1)*(NDY_BND-1)*2
          DO J = 0, NDY_BND-2
            DO I = 0, NDX_BND-2
              WRITE (3010, '(5I10)')  I+J*NDX_BND           , NDX_BND+I+J* NDX_BND, NDX_BND+I+1+J*NDX_BND, 0, 0
              WRITE (3010, '(5I10)')  NDX_BND+I+1+J*NDX_BND, I+1+J*NDX_BND        , I+J*NDX_BND          , 0, 0
            END DO
          END DO
          CLOSE(3010)
        END IF
# ifdef MPI_PARALL_GRID
      END IF
# endif
# ifdef MPI_PARALL_GRID
      IF (.NOT. MULTIPLE_IN_BOUND) THEN
        IF (myrank.eq.0) THEN
          iVect(1)=NUM_NETCDF_FILES_BND
          iVect(2)=MAX_NDT_BND_FILE
          iVect(3)=NDX_BND
          iVect(4)=NDY_BND
          DO IPROC=2,nproc
            CALL MPI_SEND(iVect,4,itype, iProc-1, 811, comm, ierr)
          END DO
          rVect(1)=OFFSET_X_BND
          rVect(2)=OFFSET_Y_BND
          rVect(3)=DX_BND
          rVect(4)=DY_BND
          DO IPROC=2,nproc
            CALL MPI_SEND(rVect,4,rtype, iProc-1, 820, comm, ierr)
          END DO
          eSize=NUM_NETCDF_FILES_BND*MAX_NDT_BND_FILE
          DO IPROC=2,nproc
            CALL MPI_SEND(NDT_BND_FILE,NUM_NETCDF_FILES_BND,itype, iProc-1, 812, comm, ierr)
            CALL MPI_SEND(BND_TIME_ALL_FILES,eSize,rtype, iProc-1, 813, comm, ierr)
          END DO
        ELSE
          CALL MPI_RECV(iVect,4,itype, 0, 811, comm, istatus, ierr)
          NUM_NETCDF_FILES_BND=iVect(1)
          MAX_NDT_BND_FILE=iVect(2)
          NDX_BND=iVect(3)
          NDY_BND=iVect(4)
          CALL MPI_RECV(rVect,4,rtype, 0, 820, comm, istatus, ierr)
          OFFSET_X_BND=rVect(1)
          OFFSET_Y_BND=rVect(2)
          DX_BND=rVect(3)
          DY_BND=rVect(4)
          ALLOCATE(BND_TIME_ALL_FILES(NUM_NETCDF_FILES_BND,MAX_NDT_BND_FILE), NDT_BND_FILE(NUM_NETCDF_FILES_BND), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 14')
          CALL MPI_RECV(NDT_BND_FILE,NUM_NETCDF_FILES_BND,itype, 0, 812, comm, istatus, ierr)
          eSize=NUM_NETCDF_FILES_BND*MAX_NDT_BND_FILE
          CALL MPI_RECV(BND_TIME_ALL_FILES,eSize,rtype, 0, 813, comm, istatus, ierr)
        END IF
      END IF
# endif
      SEBO%DELT = (BND_TIME_ALL_FILES(1,2) - BND_TIME_ALL_FILES(1,1)) * DAY2SEC
!      write(STAT%FHNDL,*) SEBO%DELT, BND_TIME_ALL_FILES(1,2), BND_TIME_ALL_FILES(1,1)
      ALLOCATE (HS_WW3(NDX_BND,NDY_BND), FP_WW3(NDX_BND,NDY_BND), T02_WW3(NDX_BND,NDY_BND), DSPR_WW3(NDX_BND,NDY_BND), DIR_WW3(NDX_BND,NDY_BND), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 15')
      HS_WW3 = 0.
      FP_WW3 = 0.
      T02_WW3 = 0.
      DSPR_WW3 = 0.
      DIR_WW3 = 0.
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef NCDF
      SUBROUTINE READ_NETCDF_WW3_IVAR(IFILE, IT, IVAR, VAR_READ)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
      integer, intent(in) :: IFILE, IT, IVAR
      character (len = *), parameter :: CallFct = "READ_NETCDF_WW3_IVAR"
      CHARACTER(LEN=40)  :: EVAR
      REAL(rkind), intent(out) :: VAR_READ(NDX_BND,NDY_BND)
      integer :: ITMP(NDX_BND,NDY_BND)
      REAL(rkind) :: scale_factor
      integer ncid, var_id
      IF (IVAR .eq. 3) THEN
        EVAR='hs'
      ELSE IF (IVAR .eq. 2) THEN
        EVAR='fp'
      ELSE IF (IVAR .eq. 5) THEN
        EVAR='t02'
      ELSE IF (IVAR .eq. 4) THEN
        EVAR='spr'
      ELSE IF (IVAR .eq. 1) THEN
        EVAR='dir'
      ELSE
        Print *, 'IVAR=', IVAR
        CALL WWM_ABORT('Wrong IVAR')
      END IF

      CALL TEST_FILE_EXIST_DIE("Missing ww3 boundary condition file : ", TRIM(NETCDF_FILE_NAMES_BND(IFILE,IVAR)))
      ISTAT = NF90_OPEN(NETCDF_FILE_NAMES_BND(IFILE,IVAR),NF90_NOWRITE,ncid)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)

      ISTAT = nf90_inq_varid(ncid, TRIM(EVAR), var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)

      ISTAT = nf90_get_att(ncid, var_id, 'scale_factor', scale_factor)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, ISTAT)

      ISTAT = NF90_GET_VAR(ncid, var_id, ITMP,  start = (/ 1, 1, IT /), count = (/ NDX_BND, NDY_BND, 1/))
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, ISTAT)
      VAR_READ = MyREAL(ITMP) * scale_factor

      ISTAT = nf90_close(ncid)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef NCDF
      SUBROUTINE READ_NETCDF_WW3_SINGLE(IFILE,IT)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IFILE, IT
      CALL READ_NETCDF_WW3_IVAR(IFILE, IT, 3, HS_WW3)
      CALL READ_NETCDF_WW3_IVAR(IFILE, IT, 2, FP_WW3)
      CALL READ_NETCDF_WW3_IVAR(IFILE, IT, 5, T02_WW3)
      CALL READ_NETCDF_WW3_IVAR(IFILE, IT, 4, DSPR_WW3)
      CALL READ_NETCDF_WW3_IVAR(IFILE, IT, 1, DIR_WW3)
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef NCDF
      SUBROUTINE READ_NETCDF_WW3_PARAM
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER     :: IFILE, IT
      INTEGER     :: counter, ip, i, j
      REAL(rkind), ALLOCATABLE    :: U(:), V(:), H(:)
      REAL(rkind), SAVE           :: TIME, scale_factor
      INTEGER IX, IY, IPROC
      REAL(rkind), ALLOCATABLE :: ARR_send_recv(:)
      integer, allocatable :: bnd_rqst(:), bnd_stat(:,:)
      CALL COMPUTE_IFILE_IT(IFILE, IT)
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_BOUND) THEN
        CALL READ_NETCDF_WW3_SINGLE(IFILE,IT)
      ELSE
        allocate(ARR_send_recv(5*NDX_BND*NDY_BND))
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 16')
        IF (myrank .eq. 0) THEN
          CALL READ_NETCDF_WW3_SINGLE(IFILE,IT)
          J=0
          DO IX=1,NDX_BND
            DO IY=1,NDY_BND
              J=J+1
              ARR_send_recv(J)=HS_WW3(IX,IY)
              J=J+1
              ARR_send_recv(J)=FP_WW3(IX,IY)
              J=J+1
              ARR_send_recv(J)=T02_WW3(IX,IY)
              J=J+1
              ARR_send_recv(J)=DSPR_WW3(IX,IY)
              J=J+1
              ARR_send_recv(J)=DIR_WW3(IX,IY)
            END DO
          END DO
          allocate(bnd_rqst(nproc-1), bnd_stat(MPI_STATUS_SIZE,nproc-1), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 17')
          DO iProc=2,nproc
            CALL mpi_isend(ARR_send_recv, 5*NDX_BND*NDY_BND, rtype, iProc-1, 2032, comm, bnd_rqst(iProc-1), ierr)
          END DO
          IF (nproc > 1) THEN
            CALL MPI_WAITALL(nproc-1, bnd_rqst, bnd_stat, ierr)
          END IF
          deallocate(bnd_rqst, bnd_stat)
        ELSE
          CALL MPI_RECV(ARR_send_recv, 5*NDX_BND*NDY_BND, rtype, 0, 2032, comm, istatus, ierr)
          J=0
          DO IX=1,NDX_BND
            DO IY=1,NDY_BND
              J=J+1
              HS_WW3(IX,IY) = ARR_send_recv(J)
              J=J+1
              FP_WW3(IX,IY) = ARR_send_recv(J)
              J=J+1
              T02_WW3(IX,IY) = ARR_send_recv(J)
              J=J+1
              DSPR_WW3(IX,IY) = ARR_send_recv(J)
              J=J+1
              DIR_WW3(IX,IY) = ARR_send_recv(J)
            END DO
          END DO
        END IF
      END IF
# else
      CALL READ_NETCDF_WW3_SINGLE(IFILE,IT)
# endif
      IF (LWRITE_WW3_RESULTS) THEN
        OPEN(3012, FILE  = 'ergwiii.bin', FORM = 'UNFORMATTED')
        ALLOCATE(U(NDX_BND*NDY_BND), V(NDX_BND*NDY_BND), H(NDX_BND*NDY_BND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 18')
        COUNTER = 1
        DO J = 1, NDY_BND
          DO I = 1, NDX_BND
            U(COUNTER) = HS_WW3(I,J)
            V(COUNTER) = DIR_WW3(I,J)
            H(COUNTER) = DSPR_WW3(I,J)
            COUNTER = COUNTER + 1
          END DO
        END DO
        TIME = TIME + 1.
        WRITE(3012) TIME
        WRITE(3012) (U(IP), V(IP), H(IP), IP = 1, NDX_BND*NDY_BND)
        DEALLOCATE(U,V,H)
      END IF
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SPPARM_INTER_STRUCT(NDX,NDY,DX,DY,OFFSET_X,OFFSET_Y, MNPT, XPT, YPT, VAL, DOPEAK)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDX, NDY
      REAL(rkind), INTENT(IN)    :: DX, DY, OFFSET_X, OFFSET_Y
      INTEGER, INTENT(IN)        :: MNPT
      REAL(rkind), INTENT(IN)    :: XPT(MNPT), YPT(MNPT)
!JZ: error
      REAL(rkind), INTENT(OUT)   :: VAL(8,MNPT) !IWBMNP)
      LOGICAL, INTENT(IN)        :: DOPEAK
      REAL(rkind) :: eVect(5)
      REAL(rkind) :: WX, WX1, WX2, WX3, WX4
      REAL(rkind) :: eVAL, HX1, HX2, LEN_X, LEN_Y
      REAL(rkind) :: DELTA_X, DELTA_Y
      INTEGER :: IVAR, I, J, IDX1, IDX2
      INTEGER :: IP, J_INT, I_INT
      DO IP = 1, MNPT
        LEN_X = XPT(IP) - OFFSET_X
        LEN_Y = YPT(IP) - OFFSET_Y
        I_INT = INT( LEN_X/DX ) + 1
        J_INT = INT( LEN_Y/DY ) + 1
        DELTA_X   = LEN_X - (I_INT - 1) * DX ! Abstand X u. Y
        DELTA_Y   = LEN_Y - (J_INT - 1) * DY !

        DO IVAR=1,5
          DO I=0,1
            DO J=0,1
              IDX1=I_INT + I
              IDX2=J_INT + J
              IF (IVAR .eq. 1) THEN
                WX=HS_WW3(IDX1,IDX2)
              ELSE IF (IVAR .eq. 2) THEN
                WX=DIR_WW3(IDX1,IDX2)
              ELSE IF (IVAR .eq. 3) THEN
                WX=FP_WW3(IDX1,IDX2)
              ELSE IF (IVAR .eq. 4) THEN
                WX=T02_WW3(IDX1,IDX2)
              ELSE IF (IVAR .eq. 5) THEN
                WX=DSPR_WW3(IDX1,IDX2)
              ELSE
                CALL WWM_ABORT('Wrong IVAR')
              END IF
              IF ((I .eq. 0).and.(J .eq. 0)) THEN
                WX1=WX
              END IF
              IF ((I .eq. 0).and.(J .eq. 1)) THEN
                WX2=WX
              END IF
              IF ((I .eq. 1).and.(J .eq. 1)) THEN
                WX3=WX
              END IF
              IF ((I .eq. 1).and.(J .eq. 0)) THEN
                WX4=WX
              END IF
            END DO
          END DO
          HX1       = WX1 + (WX4-WX1)/DX * DELTA_X
          HX2       = WX2 + (WX3-WX2)/DX * DELTA_X
          IF (WX1 .LT. 0. .OR. WX2 .LT. 0. .OR. WX3 .LT. 0. .OR. WX4 .LT. 0. ) THEN
            eVAL = 0.
          ELSE
            eVAL = HX1 + (HX2-HX1)/DY * DELTA_Y
          ENDIF
          eVect(IVAR)=eVAL
        END DO
        VAL(1,IP)=eVect(1)
        IF (DOPEAK) THEN
          IF (eVect(3) .gt. TINY(1.)) THEN
            VAL(2,IP)=ONE/eVect(3)
            VAL(5,IP)  = 2.
          END IF
        ELSE
          VAL(2,IP)=eVect(4)
          VAL(5,IP)  = -2.
        END IF
        VAL(3,IP)  = eVect(2)
        VAL(4,IP)  = eVect(5)
        VAL(6,IP)  = 1.
        VAL(7,IP)  = 0.1
        VAL(8,IP)  = 3.3

!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by user (either peak or mean)
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!     SPPARM(5): spectral shape (1-4),
!         1 - Pierson-Moskowitz
!         2 - JONSWAP
!         3 - BIN
!         4 - Gauss
!         positive peak (+) or mean frequency (-)
!     SPPARM(6): directional spreading in degree (1) or exponent (2)
!     SPPARM(7): gaussian width for the gauss spectrum 0.1
!     SPPARM(8): peak enhancement factor for the JONSWAP spectra 3.
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTER_STRUCT_BOUNDARY(NDX,NDY,DX,DY,OFFSET_X,OFFSET_Y,VAL)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDX, NDY
      REAL(rkind), INTENT(IN)    :: DX, DY, OFFSET_X, OFFSET_Y
      REAL(rkind), INTENT(OUT)   :: VAL(8,IWBMNP)
      REAL(rkind)  :: XPT(IWBMNP), YPT(IWBMNP)
      INTEGER :: IP
      DO IP=1,IWBMNP
        XPT(IP)=XP(IWBNDLC(IP))
        YPT(IP)=YP(IWBNDLC(IP))
      END DO
      CALL SPPARM_INTER_STRUCT(NDX,NDY,DX,DY,OFFSET_X,OFFSET_Y, IWBMNP, XPT, YPT, VAL, DOPEAK_BOUNDARY)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTER_STRUCT_DOMAIN(NDX,NDY,DX,DY,OFFSET_X,OFFSET_Y,VAL)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDX, NDY
      REAL(rkind), INTENT(IN)    :: DX, DY, OFFSET_X, OFFSET_Y
!JZ: error
      REAL(rkind), INTENT(OUT)   :: VAL(8,MNP) !IWBMNP)
      CALL SPPARM_INTER_STRUCT(NDX,NDY,DX,DY,OFFSET_X,OFFSET_Y, MNP, XP, YP, VAL, DOPEAK_BOUNDARY)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
!AR: check this for wave boundary nodes ....
      SUBROUTINE READWAVEPARFVCOM
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER         :: IP
# ifdef MPI_PARALL_GRID
      REAL(rkind)     :: RTMP
      INTEGER         :: IPP
# endif
      logical :: LOPEN2
!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by user (either peak or mean)
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!     SPPARM(5): spectral shape (1-4),
!                (1 - Pierson-Moskowitz,
!                 2 - JONSWAP,
!                 3 - BIN,
!                 4 - Gauss)
!                     negative peak (+) 
!                     or mean frequency (-)
!     SPPARM(6): directional spreading in degree (1) or exponent (2)
!     SPPARM(7): gaussian width for the gauss spectrum 0.1
!     SPPARM(8): peak enhancement factor for the JONSWAP spectra 3.3

      SPPARM(4,:) = 20.
      SPPARM(5,:) = 2.
      SPPARM(6,:) = 2.
      SPPARM(7,:) = 0.1
      SPPARM(8,:) = 3.3

      IF (LINHOM) THEN
         inquire(WAV%FHNDL,opened=LOPEN2)
         if(.not.LOPEN2) open(WAV%FHNDL,file=WAV%FNAME,status='old')
         !write(12,*)'WWM:',WAV%FHNDL,WAV%FNAME
         READ(WAV%FHNDL,*)
      END IF

# ifdef MPI_PARALL_GRID
      IPP = 0
      IF (LINHOM) THEN
        DO IP = 1, IWBMNPGL
          IF(ipgl(IWBNDGL(IP))%rank == myrank) THEN ! IF boundary nodes belong to local domain read values into boundary array
            IPP = IPP + 1
            READ (WAV%FHNDL, *) SPPARM(1,IPP), SPPARM(2,IPP), SPPARM(3,IPP)
          ELSE
            READ (WAV%FHNDL, *) RTMP, RTMP, RTMP ! ELSE ... throw them away ...
          ENDIF
        END DO
      ELSE
        READ (WAV%FHNDL, *) SPPARM(1,1), SPPARM(2,1), SPPARM(3,1)
        DO IP = 1, IWBMNPGL
          SPPARM(1:3,IP) = SPPARM(1:3,1)
        END DO
      END IF
# else
      IF (LINHOM) THEN
        DO IP = 1, IWBMNP
          READ (WAV%FHNDL, *) SPPARM(1,IP), SPPARM(2,IP), SPPARM(3,IP)
        END DO
      ELSE
        READ (WAV%FHNDL, *) SPPARM(1,1), SPPARM(2,1), SPPARM(3,1)
        DO IP = 1, IWBMNP
          SPPARM(1:3,IP) = SPPARM(1:3,1)
        END DO
      END IF
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READSPEC1D(LFIRST)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER              :: IP, IS
      LOGICAL, INTENT(IN)  :: LFIRST
!
!     Read Spectrum 1-D File ...
!     Second Line ... number of frequencies
!     Third  Line ... number of directions
!
      IF (LFIRST) THEN
        CALL TEST_FILE_EXIST_DIE('Missing wave file : ', TRIM(WAV%FNAME))
        OPEN(WAV%FHNDL, FILE = TRIM(WAV%FNAME), STATUS = 'OLD')   
        READ (WAV%FHNDL,*) WBMSC
        READ (WAV%FHNDL,*) WBMDC
        !WRITE(*,*) WBMSC, WBMDC
        CALL ALLOC_SPEC_BND
        DO IS = 1, WBMSC
          READ (WAV%FHNDL,*) SFRQ(IS,1)
        END DO
      END IF
      IF (LINHOM) THEN
        DO IP = 1, IWBMNP
          DO IS = 1, WBMSC
            READ (WAV%FHNDL, *) SPEG(IS,1,IP), SDIR(IS,1), SPRD(IS,1)
          END DO
        END DO
      ELSE
        DO IS = 1, WBMSC
          READ (WAV%FHNDL, *) SPEG(IS,1,1), SDIR(IS,1), SPRD(IS,1)
          !WRITE(*,*) SPEG(IS,1,1), SDIR(IS,1), SPRD(IS,1)
        END DO
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READSPEC2D
      USE DATAPOOL
      IMPLICIT NONE
      END SUBROUTINE
!**********************************************************************
!* READSPEC2D_WW3INIT1
!* Read the header of a WAVEWATCHIII binary spectral file to get 
!* dimensions of the spectral grid and number of output locations 
!* required for dynamic allocation.
!*
!* Called by GETWW3SPECTRA
!*
!* Authors: Aron Roland 
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)  
!**********************************************************************
      SUBROUTINE READSPEC2D_BINARY_WW3_INIT_SPEC
      USE DATAPOOL
      IMPLICIT NONE
      CHARACTER(LEN=30) :: GNAME
      CHARACTER(LEN=21) :: LABEL
      CALL TEST_FILE_EXIST_DIE('Missing wave file : ', TRIM(WAV%FNAME))
      OPEN(WAV%FHNDL,FILE=WAV%FNAME, STATUS='OLD',CONVERT='BIG_ENDIAN',FORM='UNFORMATTED')
      READ(WAV%FHNDL)LABEL, MSC_WW3,MDC_WW3, NP_WW3, GNAME
      WRITE(STAT%FHNDL,*)'START READSPEC2D_WW3_INIT_SPEC'
      WRITE(STAT%FHNDL,*)'LABEL, MSC_WW3,MDC_WW3, NP_WW3, GNAME'
      WRITE(STAT%FHNDL,*)LABEL, MSC_WW3,MDC_WW3, NP_WW3, GNAME
      CLOSE(WAV%FHNDL)
      WRITE(STAT%FHNDL,*)'DIRECTION NUMBER IN WW3 SPECTRUM:',MDC_WW3
      WRITE(STAT%FHNDL,*)'FREQUENCY NUMBER IN WW3 SPECTRUM:',MSC_WW3
      WRITE(STAT%FHNDL,'("+TRACE...",A)')'DONE READSPEC2D_WW3_INIT_SPEC'
      END SUBROUTINE
!**********************************************************************
!* READSPEC2D_WW3INIT2
!* Read the header of a WAVEWATCHIII binary spectral file to get
!* frequencies and directions of the spectral grid and starting time
!*
!* Called by GETWW3SPECTRA
!*
!* Authors: Aron Roland
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)
!**********************************************************************
      SUBROUTINE READSPEC2D_BINARY_WW3_INIT_TIME 
      USE DATAPOOL
      IMPLICIT NONE
      REAL                 :: SPECDMP(MSC_WW3,MDC_WW3)
      INTEGER              :: TMP, IFLAG, IP, TIME(2), IT
      INTEGER, ALLOCATABLE :: ITIME(:,:)
      CHARACTER(LEN=30)    :: GNAME
      CHARACTER(LEN=21)    :: LABEL
      CHARACTER(LEN=10)    :: PID
      CHARACTER(LEN=20)    :: CTIME1, CTIME2
      CHARACTER(LEN=15)    :: TIMESTRING

      REAL :: TMP1(MSC_WW3),TMP2(MDC_WW3) !GD: in ww3 binary file, reals 
      REAL :: TMPR1, TMPR2, TMPR3, TMPR4, TMPR5, TMPR6, TMPR7

      WRITE(STAT%FHNDL,*)'START READSPEC2D_BINARY_WW3_INIT_TIME'
   
      ALLOCATE(FQ_WW3(MSC_WW3), DR_WW3(MDC_WW3), XP_WW3(NP_WW3), YP_WW3(NP_WW3), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 19')
      FQ_WW3=ZERO
      DR_WW3=ZERO
      XP_WW3=ZERO
      YP_WW3=ZERO
 
      CALL TEST_FILE_EXIST_DIE('Missing wave file : ', TRIM(WAV%FNAME))
      OPEN(WAV%FHNDL, FILE = WAV%FNAME, STATUS = 'OLD',CONVERT='BIG_ENDIAN',FORM='UNFORMATTED')
      READ(WAV%FHNDL)LABEL,TMP,TMP,TMP,GNAME
      READ(WAV%FHNDL)TMP1
      FQ_WW3 = TMP1
      READ(WAV%FHNDL)TMP2
      DR_WW3 = TMP2
      MAXSTEP_WW3 = 0

      DO 
        READ(WAV%FHNDL,IOSTAT=IFLAG) TIME
        IF (IFLAG .GT. 0) THEN
          CALL WWM_ABORT('IFLAG incorrectly set 1')
        ELSE IF (IFLAG .LT. 0) THEN
          WRITE(STAT%FHNDL,*) 'END OF FILE REACHED AT 1, WHICH IS NICE'
          EXIT
        END IF
        DO IP = 1, NP_WW3 
          READ(WAV%FHNDL,IOSTAT=IFLAG) PID,TMPR1,TMPR2,TMPR3,TMPR4,TMPR5,TMPR6,TMPR7
!          WRITE(STAT%FHNDL,'(A10,7F15.4)') PID,TMPR1,TMPR2,TMPR3,TMPR4,TMPR5,TMPR6,TMPR7
          IF (IFLAG .GT. 0) THEN
            CALL WWM_ABORT('IFLAG incorrectly set 2')
          ELSE IF (IFLAG .LT. 0) THEN
            WRITE(STAT%FHNDL,*) 'END OF FILE REACHED AT 2, WHICH IS NOT GOOD'
            EXIT
          END IF
          READ(WAV%FHNDL,IOSTAT=IFLAG) SPECDMP(:,:)
          IF (IFLAG .GT. 0) THEN
            CALL WWM_ABORT('IFLAG incorrectly set 3')
          ELSE IF (IFLAG .LT. 0) THEN
            WRITE(STAT%FHNDL,*) 'END OF FILE REACHED AT 3, WHICH IS NOT GOOD'
            EXIT
          END IF
        ENDDO
        MAXSTEP_WW3 = MAXSTEP_WW3 + 1
        IF (MAXSTEP_WW3 == 1) TSTART_WW3 = TIME
      END DO

      WRITE(STAT%FHNDL,*) 'NUMBER OF BUOYS', NP_WW3
      WRITE(STAT%FHNDL,*) 'NUMBER OF TIME STEPS IN FILE', MAXSTEP_WW3 
 
      REWIND(WAV%FHNDL)

      ALLOCATE(ITIME(MAXSTEP_WW3,2), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 20')
      ITIME = 0

      READ(WAV%FHNDL)LABEL,TMP,TMP,TMP,GNAME
      READ(WAV%FHNDL)TMP1
      READ(WAV%FHNDL)TMP2
      DO IT = 1, MAXSTEP_WW3 
        READ(WAV%FHNDL,IOSTAT=IFLAG) ITIME(IT,:)
        DO IP = 1, NP_WW3
          READ(WAV%FHNDL,IOSTAT=IFLAG) PID,TMPR1,TMPR2,TMPR3,TMPR4,TMPR5,TMPR6,TMPR7
          READ(WAV%FHNDL,IOSTAT=IFLAG) SPECDMP(:,:)
        ENDDO
      END DO

      ALLOCATE (BND_TIME_ALL_FILES(1,MAXSTEP_WW3), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 21')
      BND_TIME_ALL_FILES = ZERO

      DO IT = 1, MAXSTEP_WW3 
        WRITE(CTIME1,*) ITIME(IT,1) 
        CALL LEADINGZERO(ITIME(IT,2),CTIME2)
        TIMESTRING = ADJUSTL(TRIM(CTIME1)//'.'//ADJUSTL(CTIME2))
        CALL CT2MJD(TIMESTRING, BND_TIME_ALL_FILES(1,IT))
!        WRITE(STAT%FHNDL,*) IT, TIMESTRING, BND_TIME_ALL_FILES(1,IT)
      END DO

      DTBOUND_WW3 = (BND_TIME_ALL_FILES(1,2)-BND_TIME_ALL_FILES(1,1))*DAY2SEC
      SEBO%DELT = DTBOUND_WW3
      SEBO%BMJD = BND_TIME_ALL_FILES(1,1) 
      SEBO%EMJD = BND_TIME_ALL_FILES(1,MAXSTEP_WW3) 

      CLOSE(WAV%FHNDL)
      WRITE(STAT%FHNDL,*)'MIN. FREQ. IN WW3 SPECTRUM:',FQ_WW3(1)
      WRITE(STAT%FHNDL,*)'MAX. FREQ. IN WW3 SPECTRUM:',FQ_WW3(MSC_WW3)
      WRITE(STAT%FHNDL,*)'NUMBER OF TIME STEPS',MAXSTEP_WW3
      WRITE(STAT%FHNDL,*)'TIME INCREMENT IN SPECTRAL FILE', DTBOUND_WW3
      WRITE(STAT%FHNDL,*)'FIRST TIME STEP IN WW3 SPECTRUM FILE:',BND_TIME_ALL_FILES(1,1)
      WRITE(STAT%FHNDL,*)'BEGING TIME, END TIME and DELT of wave boundary', SEBO%BMJD, SEBO%EMJD, SEBO%DELT
      WRITE(STAT%FHNDL,*)'BEGING TIME, END TIME and DELT of simulation', MAIN%BMJD, MAIN%EMJD, MAIN%DELT
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE READSPEC2D_WW3INIT2'
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE LEADINGZERO(INTIN,STRING)
! AR: some fortran magic with leading zero's ...
      INTEGER, INTENT(IN) :: INTIN
      CHARACTER(LEN=6), INTENT(OUT) :: STRING
      write( string, '(I6)' ) INTIN 
      string = repeat('0',6-len_trim(adjustl(string)))//adjustl(string)
      END SUBROUTINE
!**********************************************************************
!* READSPEC2D_WW3
!* Reads spectra in WAVEWATCHIII binary spectral file 
!*
!* Called by GETWW3SPECTRA
!*
!* Authors: Aron Roland
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)
!*
!* Remarks: GD: Need to be modified for time interpolation. Input
!*              arguments should include date for interpolation. If
!*              constant forcing is required, this time can be set to 
!*              zero and only the first step is read.  
!* Remakrs: AR: At this time the whole file is read but this should be 
!*              replaced by a direct access call to the binary file 
!*              Gulliaume can you do this? It is not that urgent.
!**********************************************************************
      SUBROUTINE READ_SPEC_BINARY_WW3_KERNEL(ISTEP,SPECOUT)
      USE DATAPOOL
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTEP

      INTEGER :: TMP
      INTEGER :: IP, IT, TIME(2)

      REAL(rkind), INTENT(OUT) :: SPECOUT(MSC_WW3,MDC_WW3,NP_WW3)

      REAL :: SPECOUT_SGLE(MSC_WW3,MDC_WW3)
      REAL :: FQ_WW3_SNGL(MSC_WW3), DR_WW3_SNGL(MDC_WW3)
      REAL :: XP_WW3_SGLE(NP_WW3), YP_WW3_SGLE(NP_WW3)
      REAL :: D(NP_WW3),UA(NP_WW3),UD(NP_WW3),CA(NP_WW3),CD2(NP_WW3)
      REAL :: m0,m1,m2,df
      INTEGER :: is,id

      CHARACTER(LEN=30) :: GNAME
      CHARACTER(LEN=21) :: LABEL
      CHARACTER(LEN=10) :: PID(NP_WW3)

      CALL TEST_FILE_EXIST_DIE('Missing wave file : ', TRIM(WAV%FNAME))
      OPEN(WAV%FHNDL, FILE = WAV%FNAME, STATUS = 'OLD',CONVERT='BIG_ENDIAN',FORM='UNFORMATTED')
      READ(WAV%FHNDL) LABEL,TMP,TMP,TMP,GNAME
      READ(WAV%FHNDL) FQ_WW3_SNGL
      READ(WAV%FHNDL) DR_WW3_SNGL
      IF(LBCSE) THEN ! non-stationary ...
        DO IT=1,MAXSTEP_WW3
          READ(WAV%FHNDL) TIME
          DO IP = 1, NP_WW3
            READ(WAV%FHNDL) PID(IP),YP_WW3_SGLE(IP),XP_WW3_SGLE(IP),D(IP),UA(IP),UD(IP),CA(IP),CD2(IP) ! As if XP and YP would change in time ... well i leave it as it is ... 
            YP_WW3(IP) = YP_WW3_SGLE(IP)*DEGRAD
            XP_WW3(IP) = XP_WW3_SGLE(IP)*DEGRAD
            READ(WAV%FHNDL)SPECOUT_SGLE(:,:)
!            write(*,*) sum(SPECOUT_SGLE)
            SPECOUT(:,:,IP) = SPECOUT_SGLE
            m0 = 0.; m1 = 0.; m2 = 0.
            DO ID = 1,MDC_WW3-1
              DO IS = 1,MSC_WW3-1
                DF = FQ_WW3(IS+1)-FQ_WW3(IS)
                M0 = M0+((SPECOUT_SGLE(IS+1,ID)+SPECOUT_SGLE(IS,ID))/TWO)*DF*DDIR_WW3
                M1 = M1+((FQ_WW3(IS+1)-FQ_WW3(IS))/TWO)*((SPECOUT_SGLE(IS+1,ID)+SPECOUT_SGLE(IS,ID))/TWO)*DF*DDIR_WW3
                M2 = M2+(((FQ_WW3(IS+1)-FQ_WW3(IS))/TWO)**TWO)*((SPECOUT_SGLE(IS+1,ID)+SPECOUT_SGLE(IS,ID))/TWO)*DF*DDIR_WW3
!                write(*,*) 4*sqrt(m0), m1, m2, df
              ENDDO
            ENDDO
          ENDDO ! IP
          IF (IT == ISTEP) EXIT ! Read the certain timestep indicated by ISTEP ...
        ENDDO ! IT 
      ELSE ! stationary ... 
        READ(WAV%FHNDL) TIME
        DO IP = 1, NP_WW3
          READ(WAV%FHNDL)PID(IP),YP_WW3_SGLE(IP),XP_WW3_SGLE(IP),D(IP), UA(IP),UD(IP),CA(IP),CD2(IP)
          YP_WW3(IP) = YP_WW3_SGLE(IP)*DEGRAD
          XP_WW3(IP) = XP_WW3_SGLE(IP)*DEGRAD
          READ(WAV%FHNDL)SPECOUT_SGLE(:,:)
          SPECOUT(:,:,IP) = SPECOUT_SGLE
        END DO
      ENDIF
      CLOSE(WAV%FHNDL)
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE READSPEC2D_WW3'
      END SUBROUTINE
!**********************************************************************
!* READ_SPEC_WW3 
!* (doing MPI exchanges if MULTIPLE_IN_BOUND = FALSE)
!* (otherwise, all node read the same data)
!*
!* Called by GET_BINARY_WW3_SPECTRA or GET_NC_WW3_SPECTRA
!*
!* Authors: Mathieu Dutour Sikiric
!**********************************************************************
      SUBROUTINE READ_SPEC_WW3(ISTEP,SPECOUT)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ISTEP
      REAL(rkind), INTENT(OUT) :: SPECOUT(MSC_WW3,MDC_WW3,NP_WW3)
      integer, allocatable :: send_rqst(:)
      integer, allocatable :: send_stat(:,:)
      integer siz, iProc
# ifndef MPI_PARALL_GRID
      IF(IBOUNDFORMAT == 3) THEN
        CALL READ_SPEC_BINARY_WW3_KERNEL(ISTEP,SPECOUT)
      ELSEIF(IBOUNDFORMAT == 6) THEN
        CALL READ_SPEC_NC_WW3_KERNEL(ISTEP,SPECOUT)
      END IF
# else
      IF (MULTIPLE_IN_BOUND) THEN
        IF(IBOUNDFORMAT == 3) THEN
          CALL READ_SPEC_BINARY_WW3_KERNEL(ISTEP,SPECOUT)
        ELSEIF(IBOUNDFORMAT == 6) THEN
          CALL READ_SPEC_NC_WW3_KERNEL(ISTEP,SPECOUT)
        END IF
      ELSE
        siz=MSC_WW3*MDC_WW3*NP_WW3
        IF (myrank .eq. 0) THEN
          allocate(send_rqst(nproc-1), send_stat(MPI_STATUS_SIZE,nproc-1), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 30')
          IF(IBOUNDFORMAT == 3) THEN
            CALL READ_SPEC_BINARY_WW3_KERNEL(ISTEP,SPECOUT)
          ELSEIF(IBOUNDFORMAT == 6) THEN
            CALL READ_SPEC_NC_WW3_KERNEL(ISTEP,SPECOUT)
          END IF
          DO iProc=2,nproc
            CALL mpi_isend(SPECOUT, siz, rtype, iProc-1, 2034, comm, send_rqst(iProc-1), ierr)
          END DO
          IF (nproc .gt. 1) THEN
            call mpi_waitall(nproc-1, send_rqst, send_stat,ierr)
          END IF
          deallocate(send_rqst, send_stat)
        ELSE
          CALL MPI_RECV(SPECOUT, siz, rtype, 0, 2034, comm, istatus, ierr)
        END IF
      END IF
# endif
      END SUBROUTINE
!**********************************************************************
!* GET_BINARY_WW3_SPECTRA
!* Read a WAVEWATCHIII binary spectral file and do time and space 
!* interpolation if required.
!*
!* Called by WAVE_BOUNDARY_CONDITION
!*
!* Authors: Aron Roland
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)
!**********************************************************************
      SUBROUTINE GET_BINARY_WW3_SPECTRA(WBACOUT)
      USE DATAPOOL, ONLY: NP_WW3, rkind, DR_WW3, DDIR_WW3, FQ_WW3, FRLOW, LNANINFCHK, DBG, FRHIGH
      USE DATAPOOL, ONLY: LINHOM, IWBNDLC, XP, YP, XP_WW3, YP_WW3, STAT, MSC, MDC, IWBMNP, MSC_WW3
      USE DATAPOOL, ONLY: MDC_WW3
# ifdef MPI_PARALL_GRID
      USE DATAPOOL, ONLY: XLON, YLAT
# endif
      IMPLICIT NONE
      REAL(rkind), INTENT(OUT) :: WBACOUT(MSC,MDC,IWBMNP)
      INTEGER     :: IB,IPGL,IBWW3,TIME(2),IS
      REAL(rkind) :: SPEC_WW3(MSC_WW3,MDC_WW3,NP_WW3),SPEC_WWM(MSC,MDC,NP_WW3)
      REAL(rkind) :: DIST(NP_WW3),TMP(NP_WW3), INDBWW3(NP_WW3)
      REAL(rkind) :: SPEC_WW3_TMP(MSC_WW3,MDC_WW3,NP_WW3),SPEC_WW3_UNSORT(MSC_WW3,MDC_WW3,NP_WW3)
      REAL(rkind) :: JUNK(MDC_WW3),DR_WW3_UNSORT(MDC_WW3),DR_WW3_TMP(MDC_WW3)
      REAL(rkind) :: XP_WWM,YP_WWM
      INTEGER     :: IFILE, IT
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'Begin GETWW3SPECTRA'
      CALL COMPUTE_IFILE_IT(IFILE, IT)
!
! Read spectra in file
!
      CALL READ_SPEC_WW3(IT,SPEC_WW3_UNSORT)
!
! Sort directions and carries spectra along (ww3 directions are not
! montonic)
!
      DO IBWW3 = 1, NP_WW3
        DO IS = 1,MSC_WW3
          DR_WW3_TMP=DR_WW3
          SPEC_WW3_TMP(IS,:,IBWW3) = SPEC_WW3_UNSORT(IS,:,IBWW3)
          CALL SSORT2 (DR_WW3_TMP,SPEC_WW3_TMP(IS,:,IBWW3),JUNK,MDC_WW3,2)
          SPEC_WW3(IS,:,IBWW3) = SPEC_WW3_TMP(IS,:,IBWW3)
          DR_WW3 = DR_WW3_TMP
        ENDDO
        DDIR_WW3 = DR_WW3(2) - DR_WW3(1)
!          WRITE(STAT%FHNDL,*) 'AFTER SORTING', IBWW3, SUM(SPEC_WW3(:,:,IBWW3))
      ENDDO
!
! Interpolate ww3 spectra on wwm frequency grid
! GD: at the moment 360º spanning grids are mandatory
!
      IF((FQ_WW3(1).GT.FRLOW).OR.(FQ_WW3(MSC_WW3).LT.FRHIGH)) THEN
!          WRITE(STAT%FHNDL,*)'WW3 FMIN = ',FQ_WW3(1),'WWM FMIN = ',FRLOW
!          WRITE(STAT%FHNDL,*)'WW3 FMAX = ',FQ_WW3(MSC_WW3),'WWM FMAX = ', FRHIGH
!          WRITE(STAT%FHNDL,*)'WW3 spectra does not encompass the whole WWM spectra, please carefully check if this makes sense for your simulations'
        CALL SPECTRALINT(SPEC_WW3,SPEC_WWM)
      ELSE
!          WRITE(STAT%FHNDL,*)'WW3 FMIN = ',FQ_WW3(1),'WWM FMIN = ',FRLOW
!          WRITE(STAT%FHNDL,*)'WW3 FMAX = ',FQ_WW3(MSC_WW3),'WWM FMAX = ', FRHIGH
        CALL SPECTRALINT(SPEC_WW3,SPEC_WWM)
      ENDIF

      IF (LNANINFCHK) THEN
        WRITE(DBG%FHNDL,*) 'SUMS AFTER INTERPOLATION', SUM(SPEC_WW3),SUM(SPEC_WWM)
      ENDIF
!
! Interpolate ww3 spectra on wwm boundary nodes
! GD: ww3 forcing works until here. Some more debugging is needed
!
      TMP = 0
      IF(LINHOM) THEN !nearest-neighbour interpolation
        DO IB=1,IWBMNP
          IPGL = IWBNDLC(IB)
# ifdef SCHISM
          XP_WWM=XLON(IPGL)! * RADDEG 
          YP_WWM=YLAT(IPGL)! * RADDEG 
# else
          XP_WWM = XP(IPGL)
          YP_WWM = YP(IPGL)
# endif
          IF (NP_WW3 .GT. 1) THEN
            DO IBWW3=1,NP_WW3
              DIST(IBWW3)=SQRT((XP_WWM-XP_WW3(IBWW3))**2+(YP_WWM-YP_WW3(IBWW3))**2)
              INDBWW3(IBWW3)=IBWW3
            ENDDO
            CALL SSORT2 (DIST, INDBWW3, TMP, NP_WW3, 2)
            CALL SHEPARDINT2D(2, 1./DIST(1:2),MSC,MDC,SPEC_WWM(:,:,INT(INDBWW3(1:2))), WBACOUT(:,:,IB), 1)
            !WRITE(STAT%FHNDL,'(A20, 2F20.5,3F30.10)') ' AFTER INTERPOLATION ', INDBWW3(1), INDBWW3(2), sum(SPEC_WWM(:,:,INT(INDBWW3(1)))), sum(SPEC_WWM(:,:,INT(INDBWW3(2)))), SUM(WBACOUT(:,:,IB))
          ELSE
            WBACOUT(:,:,IB) = SPEC_WWM(:,:,1)
          ENDIF
        ENDDO ! IB 
      ELSE
        DO IB=1,IWBMNP
          WBACOUT(:,:,IB) = SPEC_WWM(:,:,1) 
        ENDDO
      ENDIF
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE GETWW3SPECTRA'
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_BINARY_WW3_SPECTRA 
      USE DATAPOOL
      IMPLICIT NONE
!
! Read header to get grid dimension, frequencies and directions, 
! output locations and first time step
!
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'START WITH INIT_BINARY_WW3_SPECTRA'

      CALL READSPEC2D_BINARY_WW3_INIT_SPEC
      CALL READSPEC2D_BINARY_WW3_INIT_TIME

      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE WITH INIT_BINARY_WW3_SPECTRA'
      END SUBROUTINE
!**********************************************************************
!* INIT_NC_WW3_SPECTRA
!* Read the header of a WAVEWATCHIII netcdf spectral file to get 
!* dimensions of the spectral grid and number of output locations 
!* required for dynamic allocation.
!*
!* Called by INIT_WAVE_BOUNDARY_CONDITION
!*
!* Authors: Kevin Martins (07/2018)  
!**********************************************************************
#ifdef NCDF
      SUBROUTINE INIT_NC_WW3_SPECTRA 
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
      INTEGER :: NCID
      INTEGER :: FREQ_DIM_ID, DIRS_DIM_ID, TIME_DIM_ID, STAT_DIM_ID
      INTEGER :: FREQ_VAR_ID, DIRS_VAR_ID, TIME_VAR_ID, STAT_LAT_ID, STAT_LON_ID
      REAL(rkind) :: eJD_WW3, eJD_WWM

      !------------------------------------
      ! Reading grid and station dimensions
      !------------------------------------
      ! Opening netcdf file
      CALL TEST_FILE_EXIST_DIE("Missing WW3 boundary file : ", TRIM(WAV%FNAME))
      ISTAT = NF90_OPEN(WAV%FNAME,NF90_NOWRITE,NCID)
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE OPENING NETCDF FILE",1,ISTAT)

      ! Reading frequency dimension
      ISTAT = NF90_INQ_DIMID(NCID,'frequency',FREQ_DIM_ID)
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING FREQUENCY DIMENSION ID",2,ISTAT)
      ISTAT = NF90_INQUIRE_DIMENSION(NCID,FREQ_DIM_ID,len = MSC_WW3)
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING FREQUENCY DIMENSION",3,ISTAT)

      ! Reading direction dimension
      ISTAT = NF90_INQ_DIMID(NCID,'direction',DIRS_DIM_ID)
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING DIRECTION DIMENSION ID",4,ISTAT)
      ISTAT = NF90_INQUIRE_DIMENSION(NCID,DIRS_DIM_ID,len = MDC_WW3)
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING DIRECTION DIMENSION",5,ISTAT)

      ! Reading time dimension
      ISTAT = NF90_INQ_DIMID(NCID,'time',TIME_DIM_ID)
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING DIRECTION DIMENSION ID",6,ISTAT)
      ISTAT = NF90_INQUIRE_DIMENSION(NCID,TIME_DIM_ID,len = MAXSTEP_WW3)
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING DIRECTION DIMENSION",7,ISTAT)

      ! Reading station dimension
      ISTAT = NF90_INQ_DIMID(NCID,'station',STAT_DIM_ID)
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING STATION DIMENSION ID",8,ISTAT)
      ISTAT = NF90_INQUIRE_DIMENSION(NCID,STAT_DIM_ID,len = NP_WW3)
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING STATION DIMENSION",9,ISTAT)

      !------------------------------------
      ! Reading variables
      !------------------------------------
      ! Arrays for the WW3 spectra
      ALLOCATE(FQ_WW3(MSC_WW3),DR_WW3(MDC_WW3),BND_TIME_ALL_FILES(1,MAXSTEP_WW3), XP_WW3(NP_WW3), YP_WW3(NP_WW3),stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, WW3 freq/dirs allocate error')

      ! Reading longitude
      ISTAT = NF90_INQ_VARID(NCID,'longitude',STAT_LON_ID)
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING LONGITUDE VARIABLE",10,ISTAT)
      ISTAT = NF90_GET_VAR(NCID,STAT_LON_ID,XP_WW3,start=(/1,1/),count = (/NP_WW3,1/))
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING LONGITUDE",11, ISTAT)

      ! Reading latitude
      ISTAT = NF90_INQ_VARID(NCID,'latitude',STAT_LAT_ID)
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING LATITUDE VARIABLE",12,ISTAT)
      ISTAT = NF90_GET_VAR(NCID,STAT_LAT_ID,YP_WW3,start=(/1,1/),count = (/NP_WW3,1/))
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING LATITUDE",13, ISTAT)
	  
      ! Convert coordinates to RAD
      XP_WW3 = XP_WW3*DEGRAD
      YP_WW3 = YP_WW3*DEGRAD

      ! Reading frequencies
      ISTAT = NF90_INQ_VARID(NCID,'frequency',FREQ_VAR_ID)
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING FREQUENCY VARIABLE",14,ISTAT)
      ISTAT = NF90_GET_VAR(NCID,FREQ_VAR_ID,FQ_WW3,start=(/1/),count = (/MSC_WW3/))
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING FREQUENCIES",15, ISTAT)

      ! Reading directions
      ISTAT = NF90_INQ_VARID(NCID,'direction',DIRS_VAR_ID)
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING DIRECTION VARIABLE",16,ISTAT)
      ISTAT = NF90_GET_VAR(NCID,DIRS_VAR_ID,DR_WW3,start=(/1/),count = (/MDC_WW3/))
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING DIRECTIONS",17, ISTAT)
  
      ! Reading time
      ISTAT = NF90_INQ_VARID(NCID,'time',TIME_VAR_ID)
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING DIRECTION VARIABLE",18,ISTAT)
      ISTAT = NF90_GET_VAR(NCID,TIME_VAR_ID,BND_TIME_ALL_FILES,start=(/1/),count = (/MAXSTEP_WW3/))
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING DIRECTIONS",19,ISTAT)

      !------------------------------------
      ! Conventions (direction and time)
      !------------------------------------
      ! Convention for directions: care must be taken here as sometimes the directions in WW3 
      ! are expressed with the origin taken towards the East. This corresponds to the last version of WW3.
      IF(MAXVAL(DR_WW3)>300) THEN
        DR_WW3 = MOD(DR_WW3*PI/180,PI2) ! Simple conversion to RAD
      ENDIF
      ! Converting time from WW3 convention to WWM convention
      ! In WW3, time is given in Julian days, with reference from [1990,1,1,0,0,0]
      CALL DATE2JD(1858,11,17,0,0,0,eJD_WWM)
      CALL DATE2JD(1990,1,1,0,0,0,eJD_WW3)
      BND_TIME_ALL_FILES = BND_TIME_ALL_FILES + eJD_WW3 - eJD_WWM
      DTBOUND_WW3 = (BND_TIME_ALL_FILES(1,2)-BND_TIME_ALL_FILES(1,1))*DAY2SEC
      SEBO%DELT = DTBOUND_WW3
      SEBO%BMJD = BND_TIME_ALL_FILES(1,1) 
      SEBO%EMJD = BND_TIME_ALL_FILES(1,MAXSTEP_WW3) 

      ! Closing the netcdf file
      ISTAT = NF90_CLOSE(ncid)
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE CLOSING NETCDF FILE",8,ISTAT)
      END SUBROUTINE
#endif
!**********************************************************************
!* READ_SPEC_NC_WW3_KERNEL, Called by READ_SPEC_WW3

!* Reads spectra in WAVEWATCHIII netcdf spectral file 

!* Authors: Kevin MARTINS
!**********************************************************************
#ifdef NCDF
      SUBROUTINE READ_SPEC_NC_WW3_KERNEL(ISTEP,SPECOUT)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE

      REAL(rkind), INTENT(OUT) :: SPECOUT(MSC_WW3,MDC_WW3,NP_WW3)
      INTEGER, INTENT(IN) :: ISTEP
      REAL :: SPEC_WW3(MDC_WW3,MSC_WW3)
      INTEGER :: IP, NCID, ENERGY_VAR_ID, DD, FF

      ! Opening netcdf file
      ISTAT = NF90_OPEN(WAV%FNAME,NF90_NOWRITE,NCID)
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE OPENING NETCDF FILE",1,ISTAT)

      ! Loading the energy density variable ID
      ISTAT = NF90_INQ_VARID(NCID,'efth',ENERGY_VAR_ID)
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING SPECTRA DIMENSION",2,ISTAT)

      ! Reading spectra
      IF(LBCSE) THEN ! non-stationary ...
        DO IP = 1,NP_WW3
          ISTAT = NF90_GET_VAR(NCID,ENERGY_VAR_ID,SPEC_WW3,start=(/1,1,IP,ISTEP/),count = (/MDC_WW3,MSC_WW3,1,1/))
          CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING SPECTRA",3,ISTAT)
          ! Re-ordering spectrum (WW3 is [MDC_WW3,MSC_WW3] while WWM is [MSC_WW3,MDC_WW3])
          SPECOUT(:,:,IP) = TRANSPOSE(SPEC_WW3)
        ENDDO ! IP
      ELSE ! stationary ... 
        DO IP = 1,NP_WW3
          ISTAT = NF90_GET_VAR(NCID,ENERGY_VAR_ID,SPEC_WW3,start=(/1,1,IP,1/),count = (/MDC_WW3,MSC_WW3,1,1/))
          CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE READING SPECTRA",4,ISTAT)
          ! Re-ordering spectrum (WW3 is [MDC_WW3,MSC_WW3] while WWM is [MSC_WW3,MDC_WW3])
          SPECOUT(:,:,IP) = TRANSPOSE(SPEC_WW3)
        END DO
      ENDIF
      ISTAT = NF90_CLOSE(NCID)
      CALL GENERIC_NETCDF_ERROR_WWM("ERROR WHILE CLOSING NETCDF FILE",5,ISTAT)
      END SUBROUTINE
#endif
!**********************************************************************
!* GET_NC_WW3_SPECTRA
!* Read a WAVEWATCHIII netcdf spectral file and do time and space 
!* interpolation if required.
!*
!* Called by WAVE_BOUNDARY_CONDITION
!*
!* Authors: Kevin Martins (18/12/2012)  
!**********************************************************************
#ifdef NCDF
      SUBROUTINE GET_NC_WW3_SPECTRA(WBACOUT)
      USE DATAPOOL, ONLY: NP_WW3, rkind, DR_WW3, DDIR_WW3, FQ_WW3, FRLOW, LNANINFCHK, DBG, FRHIGH
      USE DATAPOOL, ONLY: LINHOM, IWBNDLC, XP, YP, XP_WW3, YP_WW3, STAT, MSC, MDC, IWBMNP, MSC_WW3
      USE DATAPOOL, ONLY: MDC_WW3, SPDIR, PI, PI2
      USE NETCDF
# ifdef MPI_PARALL_GRID
      USE DATAPOOL, ONLY: XLON, YLAT
# endif
      IMPLICIT NONE
      REAL(rkind), INTENT(OUT) :: WBACOUT(MSC,MDC,IWBMNP)
      INTEGER     :: IB,IPGL,IBWW3,TIME(2),IS
      REAL(rkind) :: SPEC_WW3(MSC_WW3,MDC_WW3,NP_WW3),SPEC_WWM(MSC,MDC,NP_WW3)
      REAL(rkind) :: DIST(NP_WW3),TMP(NP_WW3), INDBWW3(NP_WW3)
      REAL(rkind) :: SPEC_WW3_TMP(MSC_WW3,MDC_WW3,NP_WW3),SPEC_WW3_UNSORT(MSC_WW3,MDC_WW3,NP_WW3)
      REAL(rkind) :: JUNK(MDC_WW3),DR_WW3_UNSORT(MDC_WW3),DR_WW3_TMP(MDC_WW3)
      REAL(rkind) :: XP_WWM,YP_WWM
      INTEGER     :: IFILE, IT
	  
      ! Computing time step
      CALL COMPUTE_IFILE_IT(IFILE, IT)
 
      ! Read spectra in file
      CALL READ_SPEC_WW3(IT,SPEC_WW3_UNSORT)

      ! Sort directions and carries spectra along (ww3 directions are not monotonic)
      DO IBWW3 = 1, NP_WW3
        DO IS = 1,MSC_WW3
          DR_WW3_TMP=DR_WW3
          SPEC_WW3_TMP(IS,:,IBWW3) = SPEC_WW3_UNSORT(IS,:,IBWW3)
          CALL SSORT2 (DR_WW3_TMP,SPEC_WW3_TMP(IS,:,IBWW3),JUNK,MDC_WW3,2)
          SPEC_WW3(IS,:,IBWW3) = SPEC_WW3_TMP(IS,:,IBWW3)
          DR_WW3 = DR_WW3_TMP
        ENDDO
      ENDDO

      ! Interpolate ww3 spectra on WWM frequency grid
      ! GD: at the moment 360º spanning grids are mandatory
      IF((FQ_WW3(1).GT.FRLOW).OR.(FQ_WW3(MSC_WW3).LT.FRHIGH)) THEN
!          WRITE(STAT%FHNDL,*)'WW3 FMIN = ',FQ_WW3(1),'WWM FMIN = ',FRLOW
!          WRITE(STAT%FHNDL,*)'WW3 FMAX = ',FQ_WW3(MSC_WW3),'WWM FMAX = ', FRHIGH
!          WRITE(STAT%FHNDL,*)'WW3 spectra does not encompass the whole WWM spectra, please carefully check if this makes sense for your simulations'
        CALL SPECTRALINT(SPEC_WW3,SPEC_WWM)
      ELSE
!          WRITE(STAT%FHNDL,*)'WW3 FMIN = ',FQ_WW3(1),'WWM FMIN = ',FRLOW
!          WRITE(STAT%FHNDL,*)'WW3 FMAX = ',FQ_WW3(MSC_WW3),'WWM FMAX = ', FRHIGH
        CALL SPECTRALINT(SPEC_WW3,SPEC_WWM)
      ENDIF

      ! Interpolate ww3 spectra on wwm boundary nodes
      ! GD: ww3 forcing works until here. Some more debugging is needed (e.g. for LINHOM = F)
      TMP = 0
      IF(LINHOM) THEN !nearest-neighbour interpolation
        DO IB=1,IWBMNP
          IPGL = IWBNDLC(IB)
# ifdef SCHISM
          XP_WWM = XLON(IPGL)! * RADDEG 
          YP_WWM = YLAT(IPGL)! * RADDEG 
# else
          XP_WWM = XP(IPGL)
          YP_WWM = YP(IPGL)
# endif
          IF (NP_WW3 .GT. 1) THEN
            DO IBWW3=1,NP_WW3
              DIST(IBWW3)=SQRT((XP_WWM-XP_WW3(IBWW3))**2+(YP_WWM-YP_WW3(IBWW3))**2)
              INDBWW3(IBWW3)=IBWW3
            ENDDO
            CALL SSORT2 (DIST, INDBWW3, TMP, NP_WW3, 2)
            CALL SHEPARDINT2D(2, 1./DIST(1:2),MSC,MDC,SPEC_WWM(:,:,INT(INDBWW3(1:2))), WBACOUT(:,:,IB), 1)
            !WRITE(STAT%FHNDL,'(A20, 2F20.5,3F30.10)') ' AFTER INTERPOLATION ', INDBWW3(1), INDBWW3(2), sum(SPEC_WWM(:,:,INT(INDBWW3(1)))), sum(SPEC_WWM(:,:,INT(INDBWW3(2)))), SUM(WBACOUT(:,:,IB))
          ELSE
            WBACOUT(:,:,IB) = SPEC_WWM(:,:,1)
          ENDIF
        ENDDO ! IB 
      ELSE
        DO IB=1,IWBMNP
          WBACOUT(:,:,IB) = SPEC_WWM(:,:,1) 
        ENDDO
      ENDIF
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE GETWW3SPECTRA'
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SHEPARDINT2D(NP,WEIGHT,D1,D2,Z,ZINT,P)
      USE DATAPOOL, ONLY: rkind
      IMPLICIT NONE
      INTEGER, INTENT(IN)      :: NP,P,D1,D2
      REAL(rkind), INTENT(IN)  :: WEIGHT(NP), Z(D1,D2,NP)
      REAL(rkind), INTENT(OUT) :: ZINT(D1,D2)
      INTEGER                  :: IP
      REAL                     :: SW
      SW=0
      ZINT=0
      DO IP=1,NP
        SW=SW+WEIGHT(IP)**P
        ZINT(:,:)=ZINT(:,:)+Z(:,:,IP)*(WEIGHT(IP)**P)
      ENDDO
      ZINT=ZINT/SW
      END SUBROUTINE
!**********************************************************************
!* SPECTRALINT
!* Interpolate spectrum on WWM spectral grid.
!*
!* Called by GETWW3SPECTRA
!*
!* Authors: Aron Roland
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)
!*
!* Remarks: GD)The spectral grid needs to be define over 360º and
!*             the min (resp. max) frequencies need to be smaller
!*             (resp. higher) than in WWM frequency grid.
!**********************************************************************
      SUBROUTINE SPECTRALINT(SPEC_WW3,SPEC_WWM)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(IN)  :: SPEC_WW3(MSC_WW3,MDC_WW3,NP_WW3)
      REAL(rkind), INTENT(OUT) :: SPEC_WWM(MSC,MDC,NP_WW3)
      REAL(rkind) :: SPEC_WW3_TMP(MSC_WW3,MDC,NP_WW3)
      REAL(rkind) :: DF, M0_WW3, M1_WW3, M2_WW3, M0_WWM, M1_WWM, M2_WWM
      INTEGER     :: IP,IS,ID
      REAL(rkind) :: JACOBIAN(MSC), AM, SM
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING SPECTRALINT'
      JACOBIAN = ONE/(SPSIG*PI2)! ENERGY / HZ -> ACTION / RAD
      SPEC_WW3_TMP = ZERO
      SPEC_WWM     = ZERO
      IF (LNANINFCHK) THEN
        WRITE(DBG%FHNDL,'(A20,I10,3F30.2)') 'BEFORE INTERPOLATION', SUM(SPEC_WW3), SUM(SPEC_WW3_TMP), SUM(SPEC_WWM) 
      ENDIF
      DO IP=1,NP_WW3
        DO IS=1,MSC_WW3
          CALL INTERLIND(MDC_WW3,MDC,DR_WW3,SPDIR,SPEC_WW3(IS,:,IP),SPEC_WW3_TMP(IS,:,IP))
        ENDDO
        DO ID=1,MDC 
          CALL INTERLIN (MSC_WW3,MSC,FQ_WW3,FR,SPEC_WW3_TMP(:,ID,IP),SPEC_WWM(:,ID,IP))
        ENDDO
        M0_WW3 = ZERO; M1_WW3 = ZERO; M2_WW3 = ZERO
        DO ID = 1,MDC_WW3
          DO IS = 1,MSC_WW3-1
            DF = FQ_WW3(IS+1)-FQ_WW3(IS)
            AM = (SPEC_WW3(IS+1,ID,IP)+SPEC_WW3(IS,ID,IP))/TWO
            SM = (FQ_WW3(IS+1)+FQ_WW3(IS))/TWO
            M0_WW3 =M0_WW3+AM*DF*DDIR_WW3
            M1_WW3 =M1_WW3+AM*SM*DF*DDIR_WW3
            M2_WW3 =M2_WW3+AM*SM**2*DF*DDIR_WW3
          ENDDO
        ENDDO
        M0_WWM = ZERO; M1_WWM = ZERO; M2_WWM = ZERO 
        DO ID = 1,MDC
          DO IS = 1,MSC-1
            DF = FR(IS+1)-FR(IS)
            AM = (SPEC_WWM(IS+1,ID,IP)+SPEC_WWM(IS,ID,IP))/TWO
            SM = (FR(IS+1)+FR(IS))/TWO
            M0_WWM =M0_WWM+AM*DF*DDIR
            M1_WWM =M1_WWM+AM*SM*DF*DDIR
            M2_WWM =M2_WWM+AM*SM**2*DF*DDIR
          ENDDO
        ENDDO
        WRITE(STAT%FHNDL,*) 'POINT NUMBER', IP
        WRITE(STAT%FHNDL,'(A10,2F20.10,A10,2F20.10)') 'M1 = ',M1_WW3, M1_WWM, 'M2 = ',M2_WW3, M2_WWM
      END DO

      IF (LNANINFCHK) THEN
        WRITE(DBG%FHNDL,'(A20,I10,3F30.2)') 'AFTER INTERPOLATION', IP, SUM(SPEC_WW3), SUM(SPEC_WW3_TMP), SUM(SPEC_WWM)
      ENDIF

! Do jacobian

      DO IP = 1, NP_WW3
        DO ID = 1, MDC
          SPEC_WWM(:,ID,IP) = SPEC_WWM(:,ID,IP) * JACOBIAN(:) ! convert to wave action in sigma,theta space
        END DO              
      END DO

      WRITE(STAT%FHNDL,*)'CHECKING INTEGRATED PARAMETERS AFTER JACOBIAN'
      DO IP = 1, NP_WW3
        M0_WW3 = ZERO; M1_WW3 = ZERO; M2_WW3 = ZERO
        DO ID = 1,MDC_WW3
          DO IS = 1,MSC_WW3-1
            DF = FQ_WW3(IS+1)-FQ_WW3(IS)
            AM = (SPEC_WW3(IS+1,ID,IP)+SPEC_WW3(IS,ID,IP))/TWO
            SM = (FQ_WW3(IS+1)+FQ_WW3(IS))/TWO
            M0_WW3 =M0_WW3+AM*DF*DDIR_WW3
            M1_WW3 =M1_WW3+AM*SM*DF*DDIR_WW3
            M2_WW3 =M2_WW3+AM*SM**2*DF*DDIR_WW3
          ENDDO
        ENDDO
        M0_WWM = ZERO; M1_WWM = ZERO; M2_WWM = ZERO
        DO ID = 1,MDC
          DO IS = 1,MSC-1
            DF = SPSIG(IS+1)-SPSIG(IS)
            SM = (SPSIG(IS+1)+SPSIG(IS))/TWO
            AM = (SPEC_WWM(IS+1,ID,IP)+SPEC_WWM(IS,ID,IP))/TWO * SM
            M0_WWM =M0_WWM+AM*DF*DDIR
            M1_WWM =M1_WWM+AM*SM*DF*DDIR
            M2_WWM =M2_WWM+AM*SM**2*DF*DDIR
          ENDDO
        ENDDO
        WRITE(STAT%FHNDL,*) 'POINT NUMBER', IP
        WRITE(STAT%FHNDL,'(A10,2F20.10,A10,2F20.10)') 'M1 = ',M1_WW3, M1_WWM, 'M2 = ',M2_WW3, M2_WWM
      END DO
      IF (LNANINFCHK) THEN
        WRITE(DBG%FHNDL,'(A20,I10,3F30.2)') 'AFTER JACOBIAN', IP, SUM(SPEC_WW3), SUM(SPEC_WW3_TMP), SUM(SPEC_WWM)
      ENDIF
      END SUBROUTINE
!**********************************************************************
!* INTERLIND
!* Interpolate vector on a 2-pi periodic axis (directions in wwm grid).
!*
!* Called by SPECTRALINT
!*
!* Authors: Aron Roland
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)
!*
!* Remarks: GD)This routine should be togeher with INTERLIN either here 
!*             or in wwm_aux.F90.
!**********************************************************************
      SUBROUTINE INTERLIND (NX1, NX2, X1, X2, Y1, Y2)
      USE DATAPOOL, ONLY : THR, rkind,PI2
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX1, NX2

      REAL(rkind), INTENT(IN)    :: X1(NX1), Y1(NX1)
      REAL(rkind), INTENT(IN)    :: X2(NX2)
      REAL(rkind), INTENT(OUT)   :: Y2(NX2)

      INTEGER             :: I, J
      REAL(rkind)         :: DX1(NX1)

      DO I = 1, NX1 - 1
        DX1(I) = X1(I+1)-X1(I)
      END DO
      DX1(NX1) = X1(1)+PI2-X1(NX1)

      DO I = 1, NX2
        DO J = 1,NX1 - 1
          IF (ABS(X2(I)-X1(J)) .LT. THR) THEN
            Y2(I) = Y1(J)
            EXIT
          ELSE IF (X2(I) .GT. X1(J) .AND. X2(I) .LT. X1(J+1)) THEN
            Y2(I) = Y1(J)+(Y1(J+1)-Y1(J))/DX1(J)*(X2(I)-X1(J))
            EXIT
          ELSE IF (X2(I) .GT. X1(NX1)) THEN
            Y2(I) = Y1(NX1)+(Y1(1)-Y1(NX1))/DX1(NX1)*(X2(I)-X1(NX1))
            EXIT
          ELSE IF (X2(I) .LT. X1(1)) THEN
            Y2(I) = Y1(NX1)+(Y1(1)-Y1(NX1))/DX1(NX1)*(X2(I)+PI2-X1(NX1))
            EXIT            
          END IF
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ALLOC_SPEC_BND()
      USE DATAPOOL
      IMPLICIT NONE
      IF (LINHOM) THEN
        ALLOCATE (SFRQ(WBMSC,IWBMNP), SPEG(WBMSC,WBMDC,IWBMNP), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 22')
        ALLOCATE (SDIR(WBMSC,IWBMNP), SPRD(WBMSC,IWBMNP), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 23')
      ELSE
        ALLOCATE (SFRQ(WBMSC,1), SPEG(WBMSC,WBMDC,1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 24')
        ALLOCATE (SDIR(WBMSC,1), SPRD(WBMSC,1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 25')
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READWAVEPARWW3
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER         :: IP
# ifdef MPI_PARALL_GRID
      INTEGER         :: IPP
      REAL(rkind)     :: RTMP
# endif

!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by user (either peak or mean)
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!     SPPARM(5): spectral shape (1-4),
!                (1 - Pierson-Moskowitz,
!                 2 - JONSWAP,
!                 3 - BIN,
!                 4 - Gauss)
!                     negative peak (+) 
!                     or mean frequency (-)
!     SPPARM(6): directional spreading in degree (1) or exponent (2)
!     SPPARM(7): gaussian width for the gauss spectrum 0.1
!     SPPARM(8): peak enhancement factor for the JONSWAP spectra 3.3

      SPPARM(4,:) = 20.
      SPPARM(5,:) = 2.
      SPPARM(6,:) = 2.
      SPPARM(7,:) = 0.1
      SPPARM(8,:) = 3.3

      IF (LINHOM) THEN
        READ(WAV%FHNDL,*)
      END IF

# ifdef MPI_PARALL_GRID
      IPP = 0
      IF (LINHOM) THEN
        DO IP = 1, IWBMNPGL
          IF(ipgl(IWBNDGL(IP))%rank == myrank) THEN ! IF boundary nodes belong to local domain read values into boundary array
            IPP = IPP + 1
            READ (WAV%FHNDL, *) SPPARM(1,IPP), SPPARM(2,IPP), SPPARM(3,IPP)
          ELSE
            READ (WAV%FHNDL, *) RTMP, RTMP, RTMP ! ELSE ... throw them away ...
          ENDIF
        END DO
      ELSE
        READ (WAV%FHNDL, *) SPPARM(1,1), SPPARM(2,1), SPPARM(3,1)
        DO IP = 1, IWBMNPGL
          SPPARM(1:3,IP) = SPPARM(1:3,1)
        END DO
      END IF
# else
      IF (LINHOM) THEN
        DO IP = 1, IWBMNP
          READ (WAV%FHNDL, *) SPPARM(1,IP), SPPARM(2,IP), SPPARM(3,IP)
        END DO
      ELSE
        READ (WAV%FHNDL, *) SPPARM(1,1), SPPARM(2,1), SPPARM(3,1)
        DO IP = 1, IWBMNP
          SPPARM(1:3,IP) = SPPARM(1:3,1)
        END DO
      END IF
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef NCDF
      SUBROUTINE WRITE_NETCDF_BOUND_HEADERS_1(FILE_NAME, nbTime, np_write, nbBound, BOUC_PARAM, BOUC_SPEC)
      USE DATAPOOL
      USE NETCDF
      implicit none
      character(len=140), intent(in) :: FILE_NAME
      integer, intent(in) :: nbTime, np_write, nbBound
      logical, intent(in) :: BOUC_PARAM, BOUC_SPEC
      !
      character (len = *), parameter :: UNITS = "units"
      integer one_dims, two_dims, three_dims, fifteen_dims
      integer mnp_dims, mne_dims, nfreq_dims, ndir_dims, eight_dims
      integer iret, var_id, ncid
      integer ntime_dims, iwbmnpgl_dims
      character (len = *), parameter :: CallFct="WRITE_NETCDF_BOUND_HEADERS_1"
      iret = nf90_create(TRIM(FILE_NAME), NF90_CLOBBER, ncid)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
      iret = nf90_def_dim(ncid, 'one', 1, one_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)
      iret = nf90_def_dim(ncid, 'two', 2, two_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)
      iret = nf90_def_dim(ncid, 'three', 3, three_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
      iret = nf90_def_dim(ncid, 'fifteen', 15, fifteen_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
      iret = nf90_def_dim(ncid, 'np_total', np_write, mnp_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
      iret = nf90_def_dim(ncid, 'nfreq', MSC, nfreq_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, iret)
      iret = nf90_def_dim(ncid, 'ndir', MDC, ndir_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, iret)
      iret = nf90_def_dim(ncid, 'eight',   8, eight_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, iret)
      iret = nf90_def_dim(ncid, 'IWBMNPGL', nbBound, iwbmnpgl_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, iret)
      !
      CALL WRITE_PARAM_1(ncid, one_dims)
      !
      CALL WRITE_NETCDF_TIME_HEADER(ncid, nbTime, ntime_dims)
      !
      iret=nf90_def_var(ncid,'IOBP',NF90_INT,(/ mnp_dims /), var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'integer')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, iret)
      iret=nf90_put_att(ncid,var_id,'description','boundary status of nodes')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, iret)
      iret=nf90_put_att(ncid,var_id,'case 0','interior point')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, iret)
      iret=nf90_put_att(ncid,var_id,'case 1','island')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 15, iret)
      iret=nf90_put_att(ncid,var_id,'case 2','dirichlet condition')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 16, iret)
      iret=nf90_put_att(ncid,var_id,'case 3','neumann condition')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 17, iret)
      iret=nf90_put_att(ncid,var_id,'case 4','unknown')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 18, iret)
      !
      iret=nf90_def_var(ncid,'IWBNDGL',NF90_INT,(/ iwbmnpgl_dims /), var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 19, iret)
      iret=nf90_put_att(ncid,var_id,'description','indices of boundary nodes')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 20, iret)
      !
      IF (BOUC_NETCDF_OUT_PARAM) THEN
        iret=nf90_def_var(ncid,'SPPARM',NF90_OUTTYPE_BOUC,(/ eight_dims, iwbmnpgl_dims, ntime_dims /), var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 21, iret)
        iret=nf90_put_att(ncid,var_id,'description','Parametric boundary condition')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 22, iret)
      END IF
      !
      IF (BOUC_NETCDF_OUT_SPECTRA) THEN
        iret=nf90_def_var(ncid,'WBAC',NF90_OUTTYPE_BOUC,(/ nfreq_dims, ndir_dims,  iwbmnpgl_dims, ntime_dims /), var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 23, iret)
        iret=nf90_put_att(ncid,var_id,'description','boundary wave action')
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 24, iret)
      END IF
      iret=nf90_close(ncid)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 25, iret)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_NETCDF_BOUND_HEADERS_2(FILE_NAME, np_write, IOBPwrite, nbBound, ListBound)
      USE DATAPOOL
      USE NETCDF
      implicit none
      character(len=140), intent(in) :: FILE_NAME
      integer, intent(in) :: np_write
      integer, intent(in) :: IOBPwrite(np_write)
      integer, intent(in) :: nbBound
      integer, intent(in) :: ListBound(nbBound)
      !
      integer ncid
      integer iret, var_id
      character (len = *), parameter :: CallFct="WRITE_NETCDF_BOUND_HEADERS_2"
      !
      iret=nf90_open(TRIM(FILE_NAME), NF90_WRITE, ncid)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
      !
      CALL WRITE_PARAM_2(ncid)
      !
      iret=nf90_inq_varid(ncid, "IOBP", var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)
      iret=nf90_put_var(ncid, var_id, IOBPwrite, start=(/1/), count =(/np_write/) )
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)
      !
      iret=nf90_inq_varid(ncid, "IWBNDGL", var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
      iret=nf90_put_var(ncid, var_id, ListBound, start=(/1/), count =(/nbBound/) )
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
      iret=nf90_close(ncid)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_NETCDF_BOUNDARY_OUTPUT
      USE DATAPOOL
      IMPLICIT NONE
      LOGICAL, SAVE :: IsFirstTime = .true.
      IF (IsFirstTime) THEN
        IsFirstTime=.FALSE.
# ifdef MPI_PARALL_GRID
        CALL SETUP_BOUNDARY_SCATTER_REDUCE_ARRAY
# endif
# ifdef MPI_PARALL_GRID
        IF (myrank .eq. rank_boundary) THEN
# endif
          IF (BOUC_NETCDF_OUT_PARAM .and. (.NOT. allocated(SPPARM_GL))) THEN
            allocate(SPPARM_GL(8,IWBMNPGL), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 26')
          END IF
          IF (BOUC_NETCDF_OUT_SPECTRA .and. (.NOT. allocated(WBAC_GL))) THEN
            allocate(WBAC_GL(MSC,MDC,IWBMNPGL), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 27')
          END IF
# ifdef MPI_PARALL_GRID
        END IF
# endif
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_NETCDF_BOUNDARY
      USE DATAPOOL
      USE NETCDF
      implicit none
      logical, save :: IsInitDone = .FALSE.
      character(len =256) :: FILE_NAME, PRE_FILE_NAME
      character (len = *), parameter :: CallFct="WRITE_NETCDF_BOUNDARY"
      integer iret, ncid, irec_dim, recs_his, var_id
      integer, save ::  ifile = 1
      integer LPOS, IP
      REAL(rkind)  :: eTimeDay
      integer POSITION_BEFORE_POINT, nbTime
      LPOS=POSITION_BEFORE_POINT(OUT_BOUC % FNAME)
      IF (OUT_STATION%IDEF.gt.0) THEN
        WRITE (PRE_FILE_NAME,10) OUT_BOUC % FNAME(1:LPOS),ifile
  10    FORMAT (a,'_',i4.4)
      ELSE
        WRITE (PRE_FILE_NAME,20) OUT_BOUC % FNAME(1:LPOS)
  20    FORMAT (a)
      ENDIF
      WRITE (FILE_NAME,30) TRIM(PRE_FILE_NAME)
  30  FORMAT (a,'.nc')
      IF (IsInitDone .eqv. .FALSE.) THEN
        IsInitDone=.TRUE.
        nbTime=-1
#ifdef MPI_PARALL_GRID        
        IF (myrank == 0) THEN
#endif
          CALL WRITE_NETCDF_BOUND_HEADERS_1(FILE_NAME, nbTime, np_total, IWBMNPGL, BOUC_NETCDF_OUT_PARAM, BOUC_NETCDF_OUT_SPECTRA)
          CALL WRITE_NETCDF_BOUND_HEADERS_2(FILE_NAME, np_total, IOBPtotal, IWBMNPGL, IWBNDGL)
#ifdef MPI_PARALL_GRID        
       END IF
#endif
      END IF
      !
      ! Getting the needed global arrays 
      !
      CALL INIT_NETCDF_BOUNDARY_OUTPUT
      IF (BOUC_NETCDF_OUT_PARAM .and. LBCWA) THEN
        CALL REDUCE_BOUNDARY_ARRAY_SPPARM
      END IF
      IF (BOUC_NETCDF_OUT_SPECTRA) THEN
        CALL REDUCE_BOUNDARY_ARRAY_WBAC
      END IF
!      WRITE(STAT%FHNDL,*) 'sum(WBAC)=', sum(WBAC)
!      WRITE(STAT%FHNDL,*) 'sum(SPPARM)=', sum(SPPARM)
!      WRITE(STAT%FHNDL,*) 'IWBMNP=', IWBMNP
!      WRITE(STAT%FHNDL,*) 'IWBMNPGL=', IWBMNPGL
!      FLUSH(STAT%FHNDL)
      !
      ! Now writing the boundary data
      !
      IF (myrank .eq. rank_boundary) THEN
        eTimeDay=MAIN%TMJD
        iret=nf90_open(TRIM(FILE_NAME), NF90_WRITE, ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
        iret=nf90_inquire(ncid, unlimitedDimId = irec_dim)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
        iret=nf90_inquire_dimension(ncid, irec_dim,len = recs_his)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, iret)
        recs_his=recs_his+1
        CALL WRITE_NETCDF_TIME(ncid, recs_his, eTimeDay)
        IF (BOUC_NETCDF_OUT_PARAM .and. LBCWA) THEN
          iret=nf90_inq_varid(ncid, 'SPPARM', var_id)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, iret)
          IF (NF90_RUNTYPE == NF90_OUTTYPE_BOUC) THEN
            iret=nf90_put_var(ncid,var_id,SPPARM_GL, start=(/1,1,recs_his/), count = (/8, IWBMNPGL,1/))
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, iret)
          ELSE
            iret=nf90_put_var(ncid,var_id,SNGL(SPPARM_GL), start=(/1,1,recs_his/), count = (/8, IWBMNPGL,1/))
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, iret)
          ENDIF
        END IF
        IF (BOUC_NETCDF_OUT_SPECTRA) THEN
!          WRITE(STAT%FHNDL,*) 'sum(WBAC_GL)=', sum(WBAC_GL)
!          WRITE(STAT%FHNDL,*) 'sum(SPPARM_GL)=', sum(SPPARM_GL)
!          FLUSH(STAT%FHNDL)
          !
          iret=nf90_inq_varid(ncid, 'WBAC', var_id)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, iret)
          IF (NF90_RUNTYPE == NF90_OUTTYPE_BOUC) THEN
            iret=nf90_put_var(ncid,var_id,WBAC_GL, start=(/1,1,1,recs_his/), count = (/MSC,MDC, IWBMNPGL,1/))
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, iret)
          ELSE
            iret=nf90_put_var(ncid,var_id,SNGL(WBAC_GL), start=(/1,1,1,recs_his/), count = (/MSC,MDC, IWBMNPGL,1/))
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, iret)
          ENDIF
        END IF
        iret=nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, iret)
      END IF
      IF (OUT_BOUC % IDEF.gt.0) THEN
        IF (recs_his .eq. OUT_BOUC % IDEF) THEN
          ifile=ifile+1
          IsInitDone = .FALSE.
        ENDIF
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_BOUNDARY_WBAC_SINGLE(IFILE, IT)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
      integer, intent(in) :: IFILE, IT
      character (len = *), parameter :: CallFct="READ_NETCDF_BOUNDARY_WBAC_SINGLE"
      integer ncid, var_id
      ISTAT = NF90_OPEN(BOUC_NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, ncid)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)
      ISTAT = nf90_inq_varid(ncid, 'WBAC', var_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)
      ISTAT = NF90_GET_VAR(ncid, var_id, WBAC_GL, start=(/1,1,1,IT/), count = (/MSC,MDC, IWBMNPGL,1/))
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, ISTAT)
      ISTAT = NF90_CLOSE(ncid)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, ISTAT)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_BOUNDARY_WBAC(WBACOUT)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(OUT)   :: WBACOUT(MSC,MDC,IWBMNP)
      integer IFILE, IT
      integer IP
      CALL COMPUTE_IFILE_IT_BOUND(IFILE, IT)
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_GRID) THEN
        CALL READ_NETCDF_BOUNDARY_WBAC_SINGLE(IFILE, IT)
        DO IP=1,IWBMNP
          WBACOUT(:,:,IP)=WBAC_GL(:,:,Indexes_boundary(IP))
        END DO
      ELSE
        IF (myrank .eq. rank_boundary) THEN
          CALL READ_NETCDF_BOUNDARY_WBAC_SINGLE(IFILE, IT)
        END IF
        CALL SCATTER_BOUNDARY_ARRAY_WBAC(WBACOUT)
      END IF
# else
      CALL READ_NETCDF_BOUNDARY_WBAC_SINGLE(IFILE, IT)
      WBACOUT=WBAC_GL
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_BOUNDARY_SPPARM_SINGLE(IFILE, IT)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
      integer, intent(in) :: IFILE, IT
      character (len = *), parameter :: CallFct="READ_NETCDF_BOUNDARY_SPPARM_SINGLE"
      integer ncid, var_id
!      Print *, 'IWBMNPGL=', IWBMNPGL
!      Print *, 'IT=', IT
      ISTAT = NF90_OPEN(BOUC_NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, ncid)
!      Print *, 'step 1'
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)
!      Print *, 'step 1'
      ISTAT = nf90_inq_varid(ncid, 'SPPARM', var_id)
!      Print *, 'step 2'
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)
!      Print *, 'step 3'
!      Print *, 'allocated(SPPARM_GL)=', allocated(SPPARM_GL)
      ISTAT = NF90_GET_VAR(ncid, var_id, SPPARM_GL, start=(/1,1,IT/), count = (/8, IWBMNPGL,1/))
!      Print *, 'step 4'
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, ISTAT)
!      Print *, 'step 5'
      ISTAT = NF90_CLOSE(ncid)
!      Print *, 'step 6'
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, ISTAT)
!      Print *, 'step 7'
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_BOUNDARY_SPPARM
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER IFILE, IT
      integer IP
      CALL COMPUTE_IFILE_IT_BOUND(IFILE, IT)
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_GRID) THEN
        CALL READ_NETCDF_BOUNDARY_SPPARM_SINGLE(IFILE, IT)
        DO IP=1,IWBMNP
          SPPARM(:,IP)=SPPARM_GL(:,Indexes_boundary(IP))
        END DO
      ELSE
        IF (myrank .eq. rank_boundary) THEN
          CALL READ_NETCDF_BOUNDARY_SPPARM_SINGLE(IFILE, IT)
        END IF
        CALL SCATTER_BOUNDARY_ARRAY_SPPARM
      END IF
# else
      CALL READ_NETCDF_BOUNDARY_SPPARM_SINGLE(IFILE, IT)
      SPPARM=SPPARM_GL
# endif
      END SUBROUTINE
!**********************************************************************
!*    accepted input files for NETCDF_IN_FILE in wwminput.nml         *
!*    NETCDF_IN_FILE = 'FileIn.nc' if FileIn.nc exists                *
!*    or FileIn_0001.nc, ...., FileIn_0002.nc if they exist           *
!**********************************************************************
      SUBROUTINE INIT_NETCDF_BOUNDARY_WWM
      USE NETCDF
      USE DATAPOOL
      IMPLICIT NONE
      integer POSITION_BEFORE_POINT, LPOS
      character(len=140) FILE_NAME
      logical LFLIVE
      INTEGER iFile, jFile, iTime, nbtime_mjd
      INTEGER dimids(2), varid, fid
      real(rkind) :: ConvertToDay
      real(rkind) :: eTimeStart, eTime, dtbound
      real(rkind), allocatable :: ListTime_mjd(:)
      character (len=100) :: eStrUnitTime
      integer idx
      character (len = *), parameter :: CallFct="INIT_NETCDF_BOUNDAY_WWM"
      IF (TRIM(NETCDF_IN_FILE) .eq. "unset") THEN
        CALL WWM_ABORT('NETCDF_IN_FILE must be set for running')
      END IF
      INQUIRE( FILE = TRIM(NETCDF_IN_FILE), EXIST = LFLIVE )
      !
      ! First determination of the file names
      !
      IF (LFLIVE ) THEN
        NUMBER_BOUC_NETCDF_FILE=1
        ALLOCATE(BOUC_NETCDF_FILE_NAMES(1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 26')
        BOUC_NETCDF_FILE_NAMES(1)=TRIM(NETCDF_IN_FILE)
      ELSE
        LPOS=POSITION_BEFORE_POINT(OUT_BOUC % FNAME)
        iFile=0
        DO
          jFile=iFile+1
          WRITE (FILE_NAME,10) NETCDF_IN_FILE(1:LPOS),jFile
          INQUIRE( FILE = TRIM(FILE_NAME), EXIST = LFLIVE )
          IF (LFLIVE .eqv. .FALSE.) THEN
            EXIT
          END IF
          iFile=jFile
        END DO
        NUMBER_BOUC_NETCDF_FILE=iFile
        IF (NUMBER_BOUC_NETCDF_FILE .eq. 0) THEN
          Print *, 'Error in INIT_NETCDF_BOUNDARY_WWM'
          Print *, 'NETCDF_IN_FILE=', TRIM(NETCDF_IN_FILE)
          Print *, 'FILE_NAME=', TRIM(FILE_NAME)
          CALL WWM_ABORT('We did not find the needed boundary files')
        END IF
        ALLOCATE(BOUC_NETCDF_FILE_NAMES(NUMBER_BOUC_NETCDF_FILE), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 27')
        DO iFile=1,NUMBER_BOUC_NETCDF_FILE
          WRITE (FILE_NAME,10) NETCDF_IN_FILE(1:LPOS),jFile
          BOUC_NETCDF_FILE_NAMES(iFile)=TRIM(FILE_NAME)
        END DO
      END IF
      !
      ! next reading the times
      !
      BOUND_NB_TIME=0
      DO iFile=1,NUMBER_BOUC_NETCDF_FILE
        ISTAT = nf90_open(TRIM(BOUC_NETCDF_FILE_NAMES(iFile)), nf90_nowrite, fid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)

        ISTAT = nf90_inq_varid(fid, "ocean_time", varid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)

        ISTAT = nf90_inquire_variable(fid, varid, dimids=dimids)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, ISTAT)

        ISTAT = nf90_inquire_dimension(fid, dimids(1), len=nbtime_mjd)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, ISTAT)

        ISTAT = nf90_close(fid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)
        BOUND_NB_TIME = BOUND_NB_TIME + nbtime_mjd
      END DO

      ALLOCATE(BOUND_LIST_IFILE(BOUND_NB_TIME), BOUND_LIST_IT(BOUND_NB_TIME), BOUND_LIST_TIME(BOUND_NB_TIME), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 28')
      idx=0
      DO iFile=1,NUMBER_BOUC_NETCDF_FILE
        ISTAT = nf90_open(TRIM(BOUC_NETCDF_FILE_NAMES(iFile)), nf90_nowrite, fid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)

        ISTAT = nf90_inq_varid(fid, "ocean_time", varid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)

        ISTAT = nf90_get_att(fid, varid, "units", eStrUnitTime)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, ISTAT)
        CALL CF_EXTRACT_TIME(eStrUnitTime, ConvertToDay, eTimeStart)

        ISTAT = nf90_inquire_variable(fid, varid, dimids=dimids)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, ISTAT)

        ISTAT = nf90_inquire_dimension(fid, dimids(1), len=nbtime_mjd)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, ISTAT)
        allocate(ListTime_mjd(nbtime_mjd), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 48')

        ISTAT = nf90_get_var(fid, varid, ListTime_mjd)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, ISTAT)

        ISTAT = nf90_close(fid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)

        DO iTime=1,nbtime_mjd
          idx=idx+1
          eTime=ListTime_mjd(iTime)*ConvertToDay + eTimeStart
          BOUND_LIST_IFILE(idx) = iFile
          BOUND_LIST_IT(idx) = iTime
          BOUND_LIST_TIME(idx) = eTime
        END DO
        DEALLOCATE(ListTime_mjd)
      END DO
      DTBOUND = (BOUND_LIST_TIME(2) - BOUND_LIST_TIME(1))*DAY2SEC
      SEBO%DELT = DTBOUND
      SEBO%BMJD = BOUND_LIST_TIME(1)
      SEBO%EMJD = BOUND_LIST_TIME(BOUND_NB_TIME)
      RETURN
  10  FORMAT (a,'_',i4.4)
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_IFILE_IT_BOUND(IFILE, IT)
      USE DATAPOOL
      IMPLICIT NONE
      integer, intent(out) :: IFILE, IT
      REAL(rkind) :: DeltaTime, DeltaT
      integer iTime      
      DeltaTime=BOUND_LIST_TIME(2) - BOUND_LIST_TIME(1)
      DeltaT=(MAIN%TMJD - BOUND_LIST_TIME(1))/DeltaTime
      iTime=NINT(DeltaT) + 1
      IFILE=BOUND_LIST_IFILE(iTime)
      IT=BOUND_LIST_IT(iTime)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXPORT_BOUC_WW3_FORMAT
      USE DATAPOOL
      IMPLICIT NONE
      CHARACTER(LEN=32), PARAMETER :: IDSTRBC='WAVEWATCH III BOUNDARY DATA FILE'
      CHARACTER(LEN=10), PARAMETER :: VERBPTBC = 'III  1.03 '
      integer nbNeumann, nbDirichlet, nbUnknown
      integer IP
      LOGICAL, SAVE :: IsFirst = .TRUE.
      REAL, allocatable :: XBPI(:), YBPI(:), RDBPI(:,:)
      INTEGER, allocatable :: IPBPI(:,:)
      REAL(rkind)    :: WVK,WVCG,WVKDEP,WVN,WVC,SPSIGLOC
      real, allocatable :: ABPIO(:)
      REAL(rkind) :: eCLATS, eCG, DEPLOC, eVal
      REAL XFR, eTH, FREQ1
      INTEGER NBI, idx, IB, I, J, NSPEC_out, IK, ITH, ISP, NK, NTH, IPglob
      INTEGER TheOut
      INTEGER TIME2(2)
      nbDirichlet=0
      nbNeumann=0
      nbUnknown=0
      DO IP=1,np_total
        IF (IOBPtotal(IP) == 2) THEN
          nbDirichlet=nbDirichlet+1
        END IF
        IF (IOBPtotal(IP) == 3) THEN
          nbNeumann=nbNeumann+1
        END IF
        IF (IOBPtotal(IP) == 4) THEN
          nbUnknown=nbUnknown+1
        END IF
      END DO
      IF ((nbNeumann .gt. 0).or.(nbUnknown .gt. 0)) THEN
# ifdef MPI_PARALL_GRID
        IF (myrank == 0) THEN
# endif
          IF (IsFirst .eqv. .TRUE.) THEN
            Print *, 'nbDirichlet=', nbDirichlet
            Print *, 'nbNeumann=', nbNeumann
            Print *, 'nbUnknown=', nbUnknown
            Print *, 'Those points will be put to 0'
          END IF
# ifdef MPI_PARALL_GRID
        END IF
# endif
      END IF
      NBI = nbDirichlet
      IF (NBI .ne. IWBMNPGL) THEN
        CALL WWM_ABORT('Code inconsistency error')
      END IF
      allocate(XBPI(NBI), YBPI(NBI), IPBPI(NBI,4), RDBPI(NBI,4), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('Error allocate XBPI/YBPI')
      idx=0
      DO IP=1,NP_TOTAL
        IF (IOBPtotal(IP) == 2) THEN
          idx=idx+1
          XBPI(idx)=MySNGL(XPtotal(IP))
          YBPI(idx)=MySNGL(YPtotal(IP))
        END IF
      END DO
      DO IB=1,NBI
        IPBPI(IB,1)=IB
        IPBPI(IB,2:4)=0
      END DO
      RDBPI(:,1)   = 1
      RDBPI(:,2:4) = 0
      CALL REDUCE_BOUNDARY_ARRAY_WBAC
# ifdef MPI_PARALL_GRID
      IF (myrank == 0) THEN
# endif
         TheOut = FHNDL_EXPORT_BOUC_WW3
!         Print *, 'IsFirst=', IsFirst
         NK    = MSC
         NTH   = MDC
         XFR   = MySNGL(SFAC)
         FREQ1 = MySNGL(FR(1))
         eTH   = MySNGL(SPDIR(1))
         IF (IsFirst .eqv. .TRUE.) THEN
!           Print *, 'Before open, case 1'
           OPEN(TheOut, FILE='nest.ww3', FORM='UNFORMATTED', status='new', action='write')
           WRITE(TheOut) IDSTRBC, VERBPTBC, NK, NTH, XFR, FREQ1, eTH, NBI
           WRITE(TheOut) (XBPI(I),I=1,NBI), (YBPI(I),I=1,NBI),                   &
                         ((IPBPI(I,J),I=1,NBI),J=1,4),                           &
                         ((RDBPI(I,J),I=1,NBI),J=1,4)
         ELSE
!           Print *, 'Before open, case 2'
           OPEN(TheOut, FILE='nest.ww3', FORM='UNFORMATTED', status='old', position='append', action='write')
         END IF
!         CALL COMPUTE_TFN(TIME2)
         NSPEC_out = NK*NTH
         allocate(ABPIO(NSPEC_out), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('Error allocate ABPIO')
         WRITE(TheOut) TIME2, NBI
         DO IB=1,NBI
           IPglob=IWBNDGL(IB)
           DEPLOC = MAX(DMIN,DEPtotal(IPglob))
           IF (LSPHE) THEN
             eCLATS = COS(DEGRAD*YPtotal(IPglob))
           ELSE
             eCLATS = 1
           END IF
           DO IK=1,NK
             CALL ALL_FROM_TABLE(SPSIGLOC,DEPLOC,WVK,WVCG,WVKDEP,WVN,WVC)
             eCG = WVCG              
             DO ITH=1,NTH
               ISP = ITH + (IK-1)*NTH
               eVal= WBAC_GL(IK,ITH,IB)*eCG/eCLATS
               ABPIO(ISP) = MySNGL(eVal)
             END DO
           END DO
           WRITE(TheOut) ABPIO
         END DO
         deallocate(ABPIO)
         CLOSE(TheOut)
# ifdef MPI_PARALL_GRID
      END IF
# endif
      deallocate(XBPI, YBPI, IPBPI, RDBPI)
      IsFirst=.FALSE.
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
