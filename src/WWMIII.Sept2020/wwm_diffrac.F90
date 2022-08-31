#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DIFFRA_SIMPLE()
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER     :: IP, IS
         REAL(rkind) :: ETOT, EWKTOT, ECGTOT, EAD
         REAL(rkind) :: TMP, AUX, TRANS_X(MNP), TRANS_Y(MNP)
         REAL(rkind) :: EWK(MNP), ECG(MNP), ENG(MNP)
         REAL(rkind) :: CGK(MNP)
         REAL(rkind) :: DXENG(MNP), DYENG(MNP), DXXEN(MNP), DYYEN(MNP), DXYEN(MNP)
         REAL(rkind) :: DXCGK(MNP), DYCGK(MNP)

         EWK = 0.
         ECG = 0.
         ENG = 0.
         CGK = 0.
         DO IP = 1, MNP
           IF (DEP(IP) .GT. DMIN) THEN 
             ETOT = ZERO
             EWKTOT = ZERO
             ECGTOT = ZERO
             DO IS = 1, MSC
               EAD = SUM(AC2(IS,:,IP))*DDIR*SIGPOW(IS,2)
               ETOT = ETOT + EAD
               EWKTOT = EWKTOT + WK(IS,IP) * EAD
               ECGTOT = ECGTOT + CG(IS,IP) * EAD
             END DO
             ETOT   = FRINTF * ETOT
             EWKTOT = FRINTF * EWKTOT
             ECGTOT = FRINTF * ECGTOT
             IF (ETOT .GT. VERYSMALL) THEN
               EWK(IP) = EWKTOT / ETOT
               ECG(IP) = ECGTOT / ETOT
               ENG(IP) = SQRT(ETOT)
             ELSE
               EWK(IP) = 0.
               ECG(IP) = 0. 
               ENG(IP) = 0.
             END IF 
             IF (EWK(IP) .GT. VERYSMALL) THEN
               CGK(IP) = ECG(IP) / EWK(IP)
             ELSE
               CGK(IP) = 0. 
             END IF 
           ELSE
             EWK(IP) = 0.
             ECG(IP) = 0.
             ENG(IP) = 0.
             CGK(IP) = 0.
           END IF
         END DO

         CALL DIFFERENTIATE_XYDIR(ENG  , DXENG, DYENG)
         CALL DIFFERENTIATE_XYDIR(DXENG, DXXEN, DXYEN)
         CALL DIFFERENTIATE_XYDIR(DYENG, DXYEN, DYYEN)         
         CALL DIFFERENTIATE_XYDIR(CGK  , DXCGK, DYCGK)

         IF (LSPHE) THEN
           TRANS_X = ONE/(DEGRAD*REARTH*COS(YP*DEGRAD))
           TRANS_Y = ONE/(DEGRAD*REARTH) 
           DXENG = DXENG * TRANS_X 
           DYENG = DYENG * TRANS_Y 
           DXXEN = DXXEN * TRANS_X**2 
           DYYEN = DYYEN * TRANS_Y**2
           DXCGK = DXCGK * TRANS_X 
           DYCGK = DYCGK * TRANS_Y 
         END IF
         
         DO IP = 1, MNP
            IF (DEP(IP) .GT. DMIN .AND. ENG(IP) .GT. VERYSMALL) THEN
              TMP = ECG(IP)*EWK(IP)*ENG(IP)
              IF (TMP > VERYSMALL) THEN
                 AUX = DXCGK(IP)*DXENG(IP)+DYCGK(IP)*DYENG(IP)+CGK(IP)*(DXXEN(IP)+DYYEN(IP))
                 TMP = AUX/TMP
              ELSE
                 TMP = ZERO
              END IF
              IF (TMP < -1.0_rkind) THEN
                 DIFRM(IP) = ONE
              ELSE
                 DIFRM(IP) = SQRT(ONE+TMP)
              END IF
            ELSE
              DIFRM(IP) = ONE
            END IF
            IF (DIFRM(IP) .GT. 1.2) DIFRM(IP) = 1.2
            IF (DIFRM(IP) .LT. 0.8) DIFRM(IP) = 0.8
         END DO

         CALL DIFFERENTIATE_XYDIR(DIFRM, DIFRX, DIFRY)

         IF (LSPHE) THEN
           DIFRX = DIFRX * TRANS_X 
           DIFRY = DIFRY * TRANS_Y 
         END IF

         IF (.TRUE.) THEN
           OPEN(555, FILE  = 'ergdiffr.bin'  , FORM = 'UNFORMATTED')
           WRITE(555) SNGL(RTIME)
           WRITE(555) (SNGL(DIFRX(IP)), SNGL(DIFRY(IP)),SNGL(DIFRM(IP))-1., IP = 1, MNP)
         END IF

         !WRITE(WWMDBG%FHNDL,*) MAXVAL(DIFRM), MAXVAL(DIFRX), MAXVAL(DIFRY)
         !WRITE(WWMDBG%FHNDL,*) MINVAL(DIFRM), MINVAL(DIFRX), MINVAL(DIFRY)
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
     SUBROUTINE BOTEFCT( EWK, DFBOT )
        USE DATAPOOL
        IMPLICIT NONE  

        REAL(rkind), INTENT(IN)    :: EWK(MNP)
        REAL(rkind), INTENT(INOUT) :: DFBOT(MNP)

        REAL(rkind) :: SLPH(MNP), CURH(MNP)
        REAL(rkind) :: DXDEP(MNP) , DYDEP(MNP) 
        REAL(rkind) :: DXXDEP(MNP), DXYDEP(MNP), DYYDEP(MNP)

        REAL(rkind) :: KH, BOTFC, BOTFS

        INTEGER     :: IP

        LOGICAL     :: MY_ISNAN, MY_ISINF

!        CALL SMOOTH( -1.1, MNP, XP, YP, DEP )  

        CALL DIFFERENTIATE_XYDIR(DEP, DXDEP ,  DYDEP)
        CALL DIFFERENTIATE_XYDIR(DXDEP , DXXDEP, DXYDEP)
        CALL DIFFERENTIATE_XYDIR(DYDEP , DXYDEP, DYYDEP)

        SLPH = DXDEP**2 + DYDEP**2
        CURH = DXXDEP + DYYDEP
        DFBOT = 0. 

        DO IP = 1, MNP
          IF (EWK(IP) < VERYSMALL) CYCLE  
          KH = EWK(IP)*DEP(IP)
          IF (KH > PI) CYCLE
          DFBOT(IP) = (BOTFC(KH)*CURH(IP)+BOTFS(KH)*EWK(IP)*SLPH(IP))*G9
          IF (MY_ISINF(DFBOT(IP)) .OR. MY_ISNAN(DFBOT(IP))) THEN
            WRITE(wwmerr,*)'DFBOT is NaN', IP,KH, CURH(IP), BOTFS(KH), BOTFC(KH), EWK(IP), SLPH(IP)
            call wwm_abort(wwmerr)
          END IF
        END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      FUNCTION BOTFC(KH)
        USE DATAPOOL
        IMPLICIT NONE
        REAL(rkind) :: BOTFC

        REAL(rkind), INTENT(IN) :: KH
        REAL(rkind)             :: DKH, AUX, AUX1
        REAL(rkind)             :: COSHKH, COSH2KH, SINHKH, SINH2KH, SINH3KH
        DKH = KH

        SINHKH  = MySINH(MIN(KDMAX,DKH))
        SINH2KH = MySINH(MIN(KDMAX,2.0_rkind*DKH))
        SINH3KH = MySINH(MIN(KDMAX,3.0_rkind*DKH))
        COSHKH  = MyCOSH(MIN(KDMAX,DKH))
        COSH2KH = MyCOSH(MIN(KDMAX,2.0_rkind*DKH))        
        AUX = -4.*KH*COSHKH+SINH3KH+SINHKH+8.*(KH**2)*SINHKH
        AUX1 = 8.*COSHKH**3*(2.*KH+SINH2KH)
        BOTFC = AUX/MAX(VERYSMALL,AUX1) - KH*MyTANH(KH)/MAX(VERYSMALL,(2.*(COSHKH)**2))

        IF (BOTFC .NE. BOTFC) THEN
          WRITE(wwmerr,*)'BOTFC is NaN Aron', KH, AUX, AUX1, MyTANH(KH)
          call wwm_abort(wwmerr)
        END IF
      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      FUNCTION BOTFS(KH)
        USE DATAPOOL, ONLY : THR, VERYSMALL, RKIND, KDMAX, wwmerr
        USE DATAPOOL, ONLY : ONE, TWO
        IMPLICIT NONE
        REAL(rkind) :: BOTFS

        REAL(rkind), INTENT(IN)  :: KH
        REAL(rkind)              :: SECH, SINH2KH, SINHKH, COSH2KH
        REAL(rkind)              :: AUX, AUX1

        REAL(rkind)              :: DKH
        
        DKH = KH

        IF (MyABS(MyCOSH(DKH)) > VERYSMALL) THEN
          SECH = 1. / MyCOSH(DKH)
        ELSE
          SECH = 0. 
        END IF

        SINHKH  = MySINH(MIN(KDMAX,DKH))
        SINH2KH = MySINH(MIN(KDMAX,2.0_rkind*DKH))
        COSH2KH = MyCOSH(MIN(KDMAX,2.0_rkind*DKH))
        AUX     = SECH**2/MAX(VERYSMALL,(6.0_rkind*(TWO*DKH+SINH2KH)**3))
        AUX1    = 8.0_rkind*(DKH**4)+16.0_rkind*(DKH**3)*SINH2KH- &
              &   9.0_rkind*(SINH2KH**2*COSH2KH+12.0_rkind*DKH*   &
              &   (ONE+TWO*(SINHKH)**4)*(DKH+SINH2KH))
        BOTFS = AUX * AUX1

        IF (BOTFS .NE.  BOTFS) THEN
          WRITE(wwmerr,*)'BOTFS is NaN', BOTFS, AUX, AUX1
          call wwm_abort(wwmerr)
        ENDIF
      END FUNCTION      
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CUREFCT( MNP, DXENG, DYENG, CURT, DFCUR )     
         USE DATAPOOL, ONLY : RKIND, wwmerr
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: MNP
         REAL(rkind), INTENT(IN)    :: DXENG(MNP), DYENG(MNP), CURT(MNP,2)
         REAL(rkind), INTENT(INOUT) :: DFCUR(MNP)

         REAL(rkind)                :: AUX(MNP), AUXX(MNP), AUXY(MNP)
         REAL(rkind)                :: DXAUXX(MNP), DYAUXY(MNP)
         
         AUX(:) = CURT(:,1) * DXENG(:) + CURT(:,2) * DYENG(:)
         AUXX(:) = AUX(:) * CURT(:,1)
         AUXY(:) = AUX(:) * CURT(:,2)        
         
         CALL DIFFERENTIATE_XYDIR(AUXX, DXAUXX, AUX)
         CALL DIFFERENTIATE_XYDIR(AUXY, DYAUXY, AUX)
         
         DFCUR(:) = DXAUXX(:) + DYAUXY(:)
      END SUBROUTINE     
!**********************************************************************
!*                                                                    *
!**********************************************************************
      FUNCTION BOTFC2(KH)
        USE DATAPOOL, ONLY : VERYSMALL, VERYLARGE, RKIND, ONE, TWO, KDMAX
        IMPLICIT NONE
        REAL(rkind) :: BOTFC2
        REAL(rkind), INTENT(IN) :: KH
        REAL(rkind)             :: AUX, AUX1 
        REAL(rkind)             :: SINHKH, COSHKH, SINH2KH, SINH3KH

        IF (KH .GT. VERYLARGE) THEN
          BOTFC2 = 0.
          RETURN
        END IF

        COSHKH  = MyCOSH(MIN(KDMAX,KH))
        SINHKH  = MySINH(MIN(KDMAX,KH))
        SINH2KH = MySINH(MIN(KDMAX,2.*KH))
        SINH3KH = MySINH(MIN(KDMAX,3.*KH))

        AUX = -4.0_rkind*KH*COSHKH + SINH3KH + SINHKH + 8.0_rkind*(KH**2)*SINHKH
        AUX1 = 8.0_rkind*COSHKH**3.0_rkind*(TWO*KH + SINH2KH)
        BOTFC2 = AUX / MAX(VERYSMALL,AUX1) - KH*MyTANH(KH) / (TWO*(COSHKH)**2)

        IF (BOTFC2 .NE. BOTFC2) THEN
           WRITE(*,*) 'BOTFC2'
           WRITE(*,*) SINHKH, COSHKH, SINH2KH, SINH3KH, KH
           WRITE(*,*) AUX, AUX1, KH*MyTANH(KH), (TWO*(COSHKH)**2)
           CALL WWM_ABORT('BOTFC2')
        ENDIF
      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      FUNCTION BOTFS2(KH)
        USE DATAPOOL, ONLY : VERySMALL, VERYLARGE, RKIND, ONE, TWO, KDMAX
        IMPLICIT NONE
        REAL(rkind) :: BOTFS2
        REAL(rkind), INTENT(IN) :: KH
        REAL(rkind)             :: SECH
        REAL(rkind)             :: AUX, AUX1, SINHKH, COSHKH
        REAL(rkind)             :: SINH2KH, SINH3KH, COSH2KH

        IF (KH .GT. VERyLARGE) THEN
          BOTFS2 = 0.
          RETURN
        END IF 

        COSHKH  = MyCOSH(MIN(KDMAX,KH))
        COSH2KH = MyCOSH(MIN(KDMAX,2*KH))
        SINHKH  = MySINH(MIN(KDMAX,KH))
        SINH2KH = MySINH(MIN(KDMAX,2.*KH))
        SINH3KH = MySINH(MIN(KDMAX,3.*KH))

        SECH = ONE / MAX(VERYSMALL,COSHKH)
        AUX = SECH**2 / MAX(VERYSMALL, (6.0_rkind*(TWO*KH + SINH2KH)**3.0_rkind))
        IF (AUX .GT. VERYSMALL) THEN
          AUX1 = 8.0_rkind*(KH**4.0_rkind) + 16.0_rkind*(KH**3.0_rkind)*SINH2KH &
          &   - 9.0_rkind*(SINH2KH)**2*COSH2KH &
          &   + 12.0_rkind*KH*(ONE + 2*SINHKH**4.0_rkind)*(KH + SINH2KH)
        ELSE
          AUX1 = 0.
        END IF
        BOTFS2 = AUX * AUX1
        IF (BOTFS2 .NE. BOTFS2) THEN
          WRITE(*,*) 'BOTFS2'
          WRITE(*,*) COSHKH, COSH2KH, SINHKH, SINH2KH, SINH3KH
          WRITE(*,*) AUX, AUX1, KH, SECH, SECH**2, COSHKH
          CALL WWM_ABORT('BOTFS2')
        END IF
      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CUREFCT2( MNP, DXENG, DYENG, CURT, DFCUR )
         USE DATAPOOL, ONLY : RKIND 
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: MNP
         REAL(rkind), INTENT(IN)    :: DXENG(MNP), DYENG(MNP), CURT(MNP,2)
         REAL(rkind), INTENT(INOUT) :: DFCUR(MNP)
         REAL(rkind)                :: AUX(MNP), AUXX(MNP), AUXY(MNP)
         REAL(rkind)                :: DXAUXX(MNP), DYAUXY(MNP)

         AUX(:) = CURT(:,1) * DXENG(:) + CURT(:,2) * DYENG(:)
         AUXX(:) = AUX(:) * CURT(:,1)
         AUXY(:) = AUX(:) * CURT(:,2)

         CALL DIFFERENTIATE_XYDIR(AUXX, DXAUXX, AUX)
         CALL DIFFERENTIATE_XYDIR(AUXY, DYAUXY, AUX)

         DFCUR(:) = DXAUXX(:) + DYAUXY(:)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE BOTEFCT2( MNP, EWK, DEP, DFBOT )
        USE DATAPOOL, ONLY : VERYSMALL, VERYLARGE, RKIND
        IMPLICIT NONE
        REAL(rkind), PARAMETER     :: GRAV = 9.81
        INTEGER, INTENT(IN)        :: MNP
        REAL(rkind), INTENT(IN)    :: EWK(MNP)
        REAL(rkind), INTENT(IN)    :: DEP(MNP)
        REAL(rkind), INTENT(INOUT) :: DFBOT(MNP)
        REAL(rkind)                :: SLPH(MNP), CURH(MNP)
        REAL(rkind)                :: DXDEP(MNP), DYDEP(MNP)
        REAL(rkind)                :: DXXDEP(MNP), DXYDEP(MNP), DYYDEP(MNP)
        REAL(rkind)                :: KH
        INTEGER                    :: IP
        REAL(rkind)                :: BOTFC2, BOTFS2

        CALL DIFFERENTIATE_XYDIR(DEP(1)  , DXDEP(1) ,  DYDEP(1))
        CALL DIFFERENTIATE_XYDIR(DXDEP(1), DXXDEP(1), DXYDEP(1))
        CALL DIFFERENTIATE_XYDIR(DYDEP(1), DXYDEP(1), DYYDEP(1))

        SLPH(:) = DXDEP(:)**2 + DYDEP(:)**2
        CURH(:) = DXXDEP(:) + DYYDEP(:)

        DO IP = 1, MNP
          KH = EWK(IP)*DEP(IP)
          DFBOT(IP) = (BOTFC2(KH)*CURH(IP)+BOTFS2(KH)*EWK(IP)*SLPH(IP))*GRAV
          IF (DFBOT(IP) .NE. DFBOT(IP)) THEN
            WRITE(*,*) 'DFBOT'
            WRITE(*,*) DFBOT(IP), BOTFC2(KH), BOTFS2(KH), CURH(IP), EWK(IP), SLPH(IP)
            CALL WWM_ABORT('DFBOT')
          ENDIF
        END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DIFFRA_EXTENDED()
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP, IS, ID
         REAL(rkind) :: WVC, WVK, WVCG, WVKDEP, WVN
         REAL(rkind) :: ETOT, EWKTOT, EWCTOT, ECGTOT, EAD
         REAL(rkind) :: DFWAV
         REAL(rkind) :: AUX
         REAL(rkind) :: DELTA
         REAL(rkind) :: EWK(MNP), EWC(MNP), ECG(MNP), ENG(MNP)
         REAL(rkind) :: CCG(MNP)
         REAL(rkind) :: DXENG(MNP), DYENG(MNP), DXXEN(MNP), DYYEN(MNP), DXYEN(MNP)
         REAL(rkind) :: DXCCG(MNP), DYCCG(MNP)
         REAL(rkind) :: DFCUR(MNP)
         REAL(rkind) :: DFBOT(MNP)
         REAL(rkind) :: DFCUT
         REAL(rkind) :: ETOTC, ETOTS, DM
         REAL(rkind) :: US
         REAL(rkind) :: CAUX, CAUX2
         REAL(rkind) :: NAUX

         DFBOT(:) = ZERO
         DFCUR(:) = ZERO

         DO IP = 1, MNP
            ETOT = ZERO
            EWKTOT = ZERO
            EWCTOT = ZERO
            ECGTOT = ZERO
            IF (DEP(IP) .GT. DMIN) THEN
              DO IS = 1, MSC
                CALL ALL_FROM_TABLE(SPSIG(IS),DEP(IP),WVK,WVCG,WVKDEP,WVN,WVC)
                EAD = SUM(AC2(IS,:,IP))*DDIR*SIGPOW(IS,2)
                ETOT = ETOT + EAD
                EWKTOT = EWKTOT + WVK *EAD
                EWCTOT = EWCTOT + WVC *EAD
                ECGTOT = ECGTOT + WVCG*EAD
              END DO
              ETOT   = FRINTF * ETOT
              EWKTOT = FRINTF * EWKTOT
              EWCTOT = FRINTF * EWCTOT
              ECGTOT = FRINTF * ECGTOT
              IF (ETOT .LT. VERYSMALL) THEN
                EWK(IP) = ZERO
                EWC(IP) = ZERO
                ECG(IP) = ZERO
                ENG(IP) = ZERO
                CCG(IP) = ZERO
              ELSE
                IF (MSC > 3) THEN
                  ETOT   = ETOT   + PTAIL(2)*EAD
                  EWKTOT = EWKTOT + PTAIL(4)*EAD
                  EWCTOT = EWCTOT + PTAIL(4)*EAD
                  ECGTOT = ECGTOT + PTAIL(4)*EAD
                END IF
                EWK(IP) = EWKTOT / ETOT
                EWC(IP) = EWCTOT / ETOT
                ECG(IP) = ECGTOT / ETOT
                ENG(IP) = SQRT(MAX(ETOT,1.0E-8_rkind))
                CCG(IP) = EWC(IP) * ECG(IP)
              END IF
            ELSE
              EWK(IP) = ZERO
              EWC(IP) = ZERO
              ECG(IP) = ZERO
              ENG(IP) = ZERO
              CCG(IP) = ZERO
           END IF
         END DO

         CALL DIFFERENTIATE_XYDIR(ENG, DXENG, DYENG)
         CALL DIFFERENTIATE_XYDIR(DXENG, DXXEN, DXYEN)
         CALL DIFFERENTIATE_XYDIR(DYENG, DXYEN, DYYEN)
         CALL DIFFERENTIATE_XYDIR(CCG, DXCCG, DYCCG)

         CALL BOTEFCT2( MNP, EWK, DEP, DFBOT )

         IF (LSTCU .OR. LSECU) CALL CUREFCT2( MNP, DXENG, DYENG, CURTXY, DFCUR )

         DO IP = 1, MNP
            AUX = CCG(IP)*EWK(IP)*EWK(IP)
            IF ( AUX*ENG(IP) .GT. VERYSMALL) THEN
               DFWAV = ( DXCCG(IP)*DXENG(IP)+DYCCG(IP)*DYENG(IP)+CCG(IP)*(DXXEN(IP)+DYYEN(IP)) ) / MAX(VERySMALL,ENG(IP))
               NAUX = ECG(IP) / MAX(VERYSMALL,EWC(IP))
               IF (LSTCU .OR. LSECU) THEN
                 ETOTC = ZERO
                 ETOTS = ZERO
                 DO ID = 1, MDC
                   EAD = SUM(AC2(:,ID,IP)*SIGPOW(:,2))*FRINTF*DDIR
                   ETOTC = ETOTC + EAD*COS(SPDIR(ID))
                   ETOTS = ETOTS + EAD*SIN(SPDIR(ID))
                 END DO
                 DM = ATAN2(ETOTS,ETOTC)
                 US = CURTXY(IP,1)*COS(DM)+CURTXY(IP,2)*SIN(DM)
                 CAUX = US / EWC(IP)
                 DFCUT = (TWO/NAUX+NAUX*CAUX)*CAUX
               ELSE
                 DFCUT = ZERO
                 CAUX = ZERO
               END IF
               CAUX2 = CAUX * CAUX
               DELTA = CAUX2*(ONE+CAUX)**2-NAUX*(CAUX2-NAUX)*(ONE+(DFWAV+DFBOT(IP)+DFCUR(IP))/AUX+DFCUT)
               IF (DELTA <= ZERO) THEN
                 DIFRM(IP) = ONE
               ELSE
                 DIFRM(IP) = ONE/(CAUX2-NAUX)*(CAUX*(ONE+CAUX)-SQRT(DELTA))
               END IF
            ELSE
               DIFRM(IP) = ONE
            END IF
            !IF (DIFRM(IP) .GT. 1.2) DIFRM(IP) = 1.2
            !IF (DIFRM(IP) .LT. 0.8) DIFRM(IP) = 0.8
         END DO
         CALL DIFFERENTIATE_XYDIR(DIFRM, DIFRX, DIFRY)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DIFFRACTION
        USE DATAPOOL
        IF (LDIFR) THEN
          IF (IDIFFR == 1 ) THEN
            CALL DIFFRA_SIMPLE
          ELSE IF (IDIFFR == 2) THEN
            CALL DIFFRA_EXTENDED
          END IF
        END IF
      END SUBROUTINE COMPUTE_DIFFRACTION
!**********************************************************************
!*                                                                    *
!**********************************************************************

