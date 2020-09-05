#include "wwm_functions.h"
!     Last change:  1     9 Jun 2004    1:44 am
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PARAMENG(IP,ACLOC,SME01,SME10,KME01,KMWAM,KMWAM2,      &
     &  WLM,URSELL,UBOT,ABRBOT,TMBOT,HS,ETOT,FP,TP,CP,KPP,LPP,DM,       &
     &  DSPR,PEAKDSPR,PEAKDM)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)

         REAL(rkind), INTENT(OUT)   :: SME01, SME10
         REAL(rkind), INTENT(OUT)   :: KME01, KMWAM, KMWAM2, WLM
         REAL(rkind), INTENT(OUT)   :: URSELL
         REAL(rkind), INTENT(OUT)   :: UBOT, ABRBOT, TMBOT
         REAL(rkind), INTENT(OUT)   :: HS, ETOT, KPP, FP, CP, DM, DSPR, TP, LPP
         REAL(rkind), INTENT(OUT)   :: PEAKDSPR, PEAKDM

         INTEGER             :: ID, IS

         REAL(rkind)                :: SINHKD2(MSC), ACTOTDS(MSC), ETOTD0S(MSC)
         REAL(rkind)                :: ETOTD1S(MSC)
         REAL(rkind)                :: ETOTQKD(MSC), ETOTKSS(MSC), ETOTSKD(MSC)
         REAL(rkind)                :: ETOTKDS(MSC), ETOTQKD2(MSC)
!         REAL(rkind)                :: ETOT_DSIG(MSC)

         REAL(rkind)                :: ACTOTDSbis(MSC), ETOTD0Sbis(MSC)
         REAL(rkind)                :: ETOTSKDbis(MSC), ETOTKDSbis(MSC)

         REAL(rkind)                :: ACTOT
         REAL(rkind)                :: DKTOT, EKTOT
         REAL(rkind)                :: ETOTC4, ETOTS4, PEAKFF
         REAL(rkind)                :: ETOT1, ESUMAC, HQUOTP, HQUOT
         REAL(rkind)                :: EAD, DS, EHFR, EFTAIL, DKTOT2
         REAL(rkind)                :: UB2, AB2, CGP, ETOTF3, ETOTF4
         REAL(rkind)                :: ETOTS, ETOTC, UB2bis, AB2bis, WVN
         REAL(rkind)                :: FF, DEG, EDI, CKTAIL, CETAIL
         REAL(rkind)                :: VEC2DEG, PPTAIL, SKK, SIG2

         KMWAM   = 10.0_rkind
         KMWAM2  = 10.0_rkind
         KME01   = 10.0_rkind
         SME01   = 10.0_rkind
         SME10   = 10.0_rkind
         HS      = ZERO
         ABRBOT  = 0.001_rkind
         UBOT    = ZERO
         TMBOT   = ZERO

         ETOT = ZERO
         EFTAIL = ONE / (PTAIL(1)-ONE)
         ESUMAC = ZERO
         DO IS = 1, MSC
           DO ID = 1, MDC
             ESUMAC = ESUMAC + ACLOC(IS,ID)
           END DO
         END DO

         IF (MSC .GE. 2) THEN
            DO ID = 1, MDC
              DO IS = 2, MSC
                 DS = SPSIG(IS) - SPSIG(IS-1)
                 EAD = ONEHALF*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS*DDIR
                 ETOT = ETOT + EAD
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
           END DO
         END IF
!
! 2do ... check the influence of ETOT on the results ... 
! this is the swan type integration which i do not like since
! it is only of 1st order ...
! I better like to use the trapezoid rule or rather simpson
! type integration, this is the next step for MSC .GT. 3
!
!         ETOT_DSIG(:)  = SUM(ACLOC(:,:),DIM=2) * SIGPOW(:,2) * FDIR 
!         ETOT1         = SUM(ETOT_DSIG)
!         ETOT1 = ETOT1 + ETOT_DSIG(MSC) * PTAIL(6) / FRINTF

!         ETOT = ETOT1
         HS = MAX(VERYSMALL, 4.0_rkind * SQRT(ETOT))

         IF (ETOT .GT. THR) THEN

            SINHKD2(:) = SINH(MIN(KDMAX,WK(:,IP)*DEP(IP)))**2
            ACTOTDS(:) = SUM(ACLOC(:,:),DIM=2) * SIGPOW(:,1) * FDIR
            ETOTD0S(:) = ACTOTDS(:) * SIGPOW(:,1)
            ETOTD1S(:) = ACTOTDS(:) * SIGPOW(:,2)
            ETOTQKD(:) = ETOTD0S(:) / SQRT(WK(:,IP))
            ETOTQKD2(:)= ETOTD0S(:) * SQRT(WK(:,IP))
            ETOTKSS(:) = ETOTD0S(:) * WK(:,IP)
            ETOTSKD(:) = ETOTD0S(:) / SINHKD2(:)
            ETOTKDS(:) = ETOTSKD(:) * SIGPOW(:,2)

            ACTOT = SUM(ACTOTDS)
            ETOT1 = SUM(ETOTD1S)
            DKTOT = SUM(ETOTQKD)
            DKTOT2= SUM(ETOTQKD2)
            EKTOT = SUM(ETOTKSS)
            UB2   = SUM(ETOTSKD)
            AB2   = SUM(ETOTKDS)

            ACTOTDSbis(:) = SUM(ACLOC(:,:)/ESUMAC,DIM=2) * SIGPOW(:,1)
            ETOTD0Sbis(:) = ACTOTDSbis(:) * SIGPOW(:,1)
            ETOTSKDbis(:) = ETOTD0Sbis(:) / SINHKD2(:)
            ETOTKDSbis(:) = ETOTSKDbis(:) * SIGPOW(:,2)
            UB2bis   = SUM(ETOTSKDbis)
            AB2bis   = SUM(ETOTKDSbis)

            ACTOT   = ACTOT  + PTAIL(5)  * ACTOTDS(MSC) / FRINTF
            ETOT1   = ETOT1  + PTAIL(7)  * ETOTD0S(MSC) * SIGPOW(MSC,1) / FRINTF
            DKTOT   = DKTOT  + PTAIL(5)  * ETOTD0S(MSC) / (SQRT(WK(MSC,IP)) * FRINTF)
            DKTOT2  = DKTOT2 + PTAIL(5)  * ETOTD0S(MSC) * (SQRT(WK(MSC,IP)) * FRINTF)
            EKTOT   = EKTOT  + PTAIL(8)  * ETOTD0S(MSC) * WK(MSC,IP) / FRINTF

            IF (ETOT > VERYSMALL) SME01  = ETOT1 / ETOT
            IF (ETOT > VERySMALL) KME01  = EKTOT / ETOT
            IF (ACTOT > VERySMALL) SME10  = ETOT / ACTOT
            IF (DKTOT > VERySMALL) KMWAM    = (ETOT/DKTOT)**2
            IF (DKTOT2 > VERySMALL) KMWAM2  = (DKTOT2/ETOT)**2
            IF (UB2   > VERYSMALL) UBOT   = SQRT(UB2)
            !IF (UB2   > SMALL) ORBITAL(IP)   = SQRT(UB2)
            IF (AB2   > VERYSMALL) ABRBOT = SQRT(2*AB2)
            IF (UB2bis .GT. THR) THEN
              IF (AB2bis/UB2bis > THR) THEN
                 TMBOT = PI2*SQRT(AB2bis/UB2bis)
              END IF
            END IF
            URSELL = (G9*HS) / (TWO*SQRT(TWO)*SME01**2*DEP(IP)**2)
         ELSE

            HS           = THR
            ABRBOT       = THR
            UBOT         = THR
            TMBOT        = THR
            SME01        = 10.0_rkind
            SME10        = 10.0_rkind
            KME01        = 10.0_rkind
            KMWAM        = 10.0_rkind
            KMWAM2       = 10.0_rkind
            URSELL       = THR 

         END IF
!
! Peak period continues version... Taken from Thesis Alves ... correct citation is given there ... :)
!
         IF (ESUMAC.gt.VERYSMALL) THEN
            ETOTF3 = ZERO
            ETOTF4 = ZERO
            ETOTC4 = ZERO
            ETOTS4 = ZERO
            DO IS = 1, MSC
               DO ID = 1, MDC
                  HQUOT=ACLOC(IS,ID)/ESUMAC
                  HQUOTP=HQUOT**4
                  ETOTF3 = ETOTF3 + SPSIG(IS) * HQUOTP * DS_BAND(IS)
                  ETOTF4 = ETOTF4 +             HQUOTP * DS_BAND(IS)
                  ETOTC4 = ETOTC4 + COSTH(ID) * HQUOTP * DS_BAND(IS)
                  ETOTS4 = ETOTS4 + SINTH(ID) * HQUOTP * DS_BAND(IS)
               END DO
            END DO
            IF(ETOTF4 .GT. VERYSMALL) THEN
               FP = ETOTF3/ETOTF4
               CALL WAVEKCG(DEP(IP), FP, WVN, CP, KPP, CGP)
               TP = ONE/(FP/PI2)
               LPP = ONE/KPP*PI2
            ELSE
               TP  = ZERO 
               FP  = ZERO
               CP  = ZERO
               KPP = 10.0_rkind
               CGP = ZERO
               LPP = ZERO
            END IF
            IF (ETOTF4 .gt. THR) THEN
               PEAKDM    = VEC2DEG (ETOTC4, ETOTS4)
               CALL DEG2NAUT(PEAKDM,DEG,LNAUTOUT)
               PEAKDM = DEG
               PEAKFF = MIN (ONE, SQRT(ETOTC4*ETOTC4+ETOTS4*ETOTS4)/ETOTF4)
               PEAKDSPR = SQRT(2.0_rkind - 2.0_rkind*PEAKFF) * 180.0_rkind/PI
            ELSE
               FF = ZERO
               PEAKDSPR = ZERO
               PEAKDM = ZERO
            END IF
         ELSE
            TP  = ZERO 
            FP  = ZERO
            CP  = ZERO
            KPP = 10.0_rkind
            CGP = ZERO
            LPP = ZERO
         END IF

!         WRITE(*,*) FP, ETOTF3, ETOTF4
         ETOTC = ZERO
         ETOTS = ZERO
         ETOT1  = ZERO
         DO ID = 1, MDC
           EAD = ZERO
           IF (MSC .GE. 2) THEN
             DO  IS = 2, MSC 
               DS  = SPSIG(IS)-SPSIG(IS-1)
               EDI = ONEHALF*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS
               EAD = EAD + EDI
            END DO
             IF (MSC .GT. 3) THEN
               EHFR = ACLOC(MSC,ID) * SPSIG(MSC)
               EAD = EAD + EHFR * SPSIG(MSC) * EFTAIL
             ENDIF
             EAD = EAD * DDIR
             ETOT1 = ETOT1 + EAD
             ETOTC  = ETOTC + EAD * COSTH(ID)
             ETOTS  = ETOTS + EAD * SINTH(ID)
           ELSE
             DS = SGHIGH - SGLOW
             EAD = ACLOC(1,ID) * DS * DDIR
             EAD = EAD * DDIR
             ETOT1 = ETOT1 + EAD
             ETOTC  = ETOTC + EAD * COSTH(ID)
             ETOTS  = ETOTS + EAD * SINTH(ID)
           END IF
         END DO

         IF (ETOT > THR ) THEN
           DM    = VEC2DEG (ETOTC, ETOTS)
           CALL DEG2NAUT(DM,DEG,LNAUTOUT)
           DM = DEG
           FF = MIN (ONE, SQRT(ETOTC*ETOTC+ETOTS*ETOTS)/ETOT)
           DSPR = SQRT(TWO-TWO*FF) * 180.0_rkind/PI
         ELSE
           FF = ZERO
           DM = ZERO
           DSPR = ZERO
         END IF

         ETOT1  = ZERO
         EKTOT = ZERO
         DO IS=1, MSC
            SIG2 = SIGPOW(IS,2)
            SKK  = SIG2 * (WK(IS,IP))**ONE!OUTPAR(3)
            DO ID=1,MDC
              ETOT1  = ETOT1 + SIG2 * ACLOC(IS,ID)
              EKTOT = EKTOT + SKK * ACLOC(IS,ID)
            ENDDO
         ENDDO
         ETOT1  = FRINTF * ETOT1
         EKTOT = FRINTF * EKTOT
         IF (MSC .GT. 3) THEN
            PPTAIL = PTAIL(1) - ONE
            CETAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))
            PPTAIL = PTAIL(1) - ONE - 2.0_rkind*ONE!OUTPAR(3)
            CKTAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))
            DO ID=1,MDC
              ETOT1   = ETOT1 + CETAIL * SIG2 * ACLOC(MSC,ID)
              EKTOT  = EKTOT + CKTAIL * SKK * ACLOC(MSC,ID)
            ENDDO
         ENDIF
         IF (ETOT.GT.ZERO) THEN
            WLM = PI2 * (ETOT1 / EKTOT) ** ONE!(1./OUTPAR(3))     
         ELSE
            WLM = ZERO
         ENDIF
      END SUBROUTINE
!**********************************************************************:
!*                                                                    *
!**********************************************************************
      SUBROUTINE STOKES_DRIFT_SURFACE_BAROTROPIC(IP,STOKESBOTTX,STOKESBOTTY,STOKESSURFX,STOKESSURFY,STOKESBAROX,STOKESBAROY)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP

         REAL(rkind), INTENT(OUT)   :: STOKESBOTTX, STOKESBOTTY
         REAL(rkind), INTENT(OUT)   :: STOKESSURFX, STOKESSURFY
         REAL(rkind), INTENT(OUT)   :: STOKESBAROX, STOKESBAROY
         INTEGER                    :: ID, IS
         REAL(rkind)                :: eQuot1, eQuot2, eProd1, eProd2
         REAL(rkind)                :: eMult, eWk, kD, eWkReal
         REAL(rkind)                :: eSinh2kd, eSinhkd, eSinhkd2
         REAL(rkind)                :: eSigma, eLoc, eUint, eVint
         REAL(rkind)                :: eDep
         REAL(rkind)                :: eProd0, eQuot0

         STOKESBOTTX=0
         STOKESBOTTY=0
         STOKESSURFX=0
         STOKESSURFY=0
         STOKESBAROX=0
         STOKESBAROY=0
         eDep=DEP(IP)
         DO IS=1,MSC
           eMult=SPSIG(IS)*DDIR*DS_INCR(IS)
           eWk=WK(IS,IP)
           kD=MIN(KDMAX, eWk*eDep)
           eWkReal=kD/eDep
           eSinh2kd=MySINH(2*kD)
           eSinhkd=MySINH(kD)
           eSinhkd2=eSinhkd**2
           eSigma=SPSIG(IS)
           eUint=0
           eVint=0
           DO ID=1,MDC
             eLoc=AC2(IS,ID,IP)*eMult
             eUint=eUint + eLoc*COSTH(ID)
             eVint=eVint + eLoc*SINTH(ID)
           END DO
           eQuot0=ONE/eSinhkd2
           eProd0=eSigma*eWkReal*eQuot0
           eQuot1=MyCOSH(2*kD)/eSinhkd2
           eProd1=eSigma*eWkReal*eQuot1
           eQuot2=(eSinh2kd/(2*kD))/eSinhkd2
           eProd2=eSigma*eWkReal*eQuot2
           STOKESBOTTX = STOKESBOTTX + eUint*eProd0
           STOKESBOTTY = STOKESBOTTY + eVint*eProd0
           STOKESSURFX = STOKESSURFX + eUint*eProd1
           STOKESSURFY = STOKESSURFY + eVint*eProd1
           STOKESBAROX = STOKESBAROX + eUint*eProd2
           STOKESBAROY = STOKESBAROY + eVint*eProd2
         END DO
      END SUBROUTINE
!**********************************************************************:
!*                                                                    *
!**********************************************************************
      SUBROUTINE STOKES_DRIFT_SURFACE_BAROTROPIC_LOC(ACLOC,DEPLOC,WKLOC,  &
     &    STOKESBOTTX,STOKESBOTTY,STOKESSURFX,STOKESSURFY,STOKESBAROX,STOKESBAROY)

         USE DATAPOOL
         IMPLICIT NONE

         REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(IN)    :: DEPLOC
         REAL(rkind), INTENT(IN)    :: WKLOC(MSC)
         REAL(rkind), INTENT(OUT)   :: STOKESBOTTX, STOKESBOTTY
         REAL(rkind), INTENT(OUT)   :: STOKESSURFX, STOKESSURFY
         REAL(rkind), INTENT(OUT)   :: STOKESBAROX, STOKESBAROY
         INTEGER             :: ID, IS
         REAL(rkind)              :: eQuot1, eQuot2, eProd1, eProd2
         REAL(rkind)              :: eMult, eWk, kD, eWkReal
         REAL(rkind)              :: eSinh2kd, eSinhkd, eSinhkd2
         REAL(rkind)              :: eSigma, eLoc, eUint, eVint
         REAL(rkind)              :: eDep
         REAL(rkind)                :: eProd0, eQuot0

         STOKESBOTTX=0
         STOKESBOTTY=0
         STOKESSURFX=0
         STOKESSURFY=0
         STOKESBAROX=0
         STOKESBAROY=0
         eDep=DEPLOC
         DO IS=1,MSC
           eMult=SPSIG(IS)*DDIR*DS_INCR(IS)
           eWk=WKLOC(IS)
           kD=MIN(KDMAX, eWk*eDep)
           eWkReal=kD/eDep
           eSinh2kd=MySINH(2*kD)
           eSinhkd=MySINH(kD)
           eSinhkd2=eSinhkd**2
           eSigma=SPSIG(IS)
           eUint=0
           eVint=0
           DO ID=1,MDC
             eLoc=ACLOC(IS,ID)*eMult
             eUint=eUint + eLoc*COSTH(ID)
             eVint=eVint + eLoc*SINTH(ID)
           END DO
           eQuot0=ONE/eSinhkd2
           eProd0=eSigma*eWkReal*eQuot0
           eQuot1=MyCOSH(2*kD)/eSinhkd2
           eProd1=eSigma*eWkReal*eQuot1
           eQuot2=(eSinh2kd/(2*kD))/eSinhkd2
           eProd2=eSigma*eWkReal*eQuot2
           STOKESBOTTX = STOKESBOTTX + eUint*eProd0
           STOKESBOTTY = STOKESBOTTY + eVint*eProd0
           STOKESSURFX = STOKESSURFX + eUint*eProd1
           STOKESSURFY = STOKESSURFY + eVint*eProd1
           STOKESBAROX = STOKESBAROX + eUint*eProd2
           STOKESBAROY = STOKESBAROY + eVint*eProd2
         END DO
      END SUBROUTINE
!**********************************************************************:
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP

         REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: SME01, SME10, KME01
         REAL(rkind), INTENT(OUT)   :: KMWAM, KMWAM2
         REAL(rkind), INTENT(OUT)   :: HS

         INTEGER             :: ID, IS

         REAL(rkind)                :: ACTOT, ETOT
         REAL(rkind)                :: ETOT_SPSIG
         REAL(rkind)                :: ETOT_WK
         REAL(rkind)                :: ETOT_ISQ_WK
         REAL(rkind)                :: ETOT_SQ_WK

         REAL(rkind)                :: Y(MSC), tmp(msc)
         REAL(rkind)                :: DS, ATAIL, ETAIL, ESIGTAIL
!         REAL(rkind)                :: dintspec, dintspec_y
!
! total energy ...
! 2do improve efficiency ... 
! 2do check integration style ... 
!
         !ETOT = DINTSPEC(IP,ACLOC)
         ETOT = ZERO
         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig 
           ETOT = ETOT + tmp(1) * ONEHALF * ds_incr(1)*ddir
           do is = 2, msc
             ETOT = ETOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
           end do
           ETOT = ETOT + ONEHALF * tmp(msc) * ds_incr(msc)*ddir
         end do
!
! if etot too small skip ...
!
         if (etot .gt. thr) then
!
! integrals ... inlined ... for speed ...
!
           ACTOT = ZERO
           y = ONE/SPSIG
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ACTOT  = ACTOT + tmp(1) * ONEHALF * ds_incr(1)*ddir
             do is = 2, msc
               ACTOT = ACTOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ACTOT  = ACTOT + ONEHALF * tmp(msc) * ds_incr(msc)*ddir
           end do
           !tmp = ONE/SPSIG
           !ACTOT       = DINTSPEC_Y(IP,ACLOC,tmp)
           ETOT_SPSIG = ZERO
           y = SIGPOW(:,1) 
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ETOT_SPSIG = ETOT_SPSIG + tmp(1) * ONEHALF * ds_incr(1)*ddir
             do is = 2, msc
               ETOT_SPSIG = ETOT_SPSIG + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ETOT_SPSIG = ETOT_SPSIG + ONEHALF * tmp(msc) * ds_incr(msc)*ddir
           end do
           !tmp = SIGPOW(:,1)
           !ETOT_SPSIG  = DINTSPEC_Y(IP,ACLOC,tmp)
           ETOT_WK = ZERO
           y = WK(:,IP) 
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ETOT_WK = ETOT_WK + ONEHALF * tmp(1) * ds_incr(1)*ddir
             do is = 2, msc
               ETOT_WK = ETOT_WK + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ETOT_WK = ETOT_WK + ONEHALF * tmp(msc) * ds_incr(msc)*ddir
           end do
           !tmp = WK(:,IP)
           !ETOT_WK     = DINTSPEC_Y(IP,ACLOC,tmp)
           ETOT_ISQ_WK = ZERO 
           y = ONE/SQRT(WK(:,IP))
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ETOT_ISQ_WK = ETOT_ISQ_WK + ONEHALF * tmp(1) * ds_incr(1)*ddir
             do is = 2, msc
               ETOT_ISQ_WK = ETOT_ISQ_WK+ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ETOT_ISQ_WK = ETOT_ISQ_WK+ONEHALF*tmp(msc) * ds_incr(msc)*ddir
           end do
           !tmp = ONE/SQRT(WK(:,IP))
           !ETOT_ISQ_WK = DINTSPEC_Y(IP,ACLOC,tmp)
           ETOT_SQ_WK = ZERO
           y = SQRT(WK(:,IP)) 
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ETOT_SQ_WK = ETOT_SQ_WK + ONEHALF * tmp(1) * ds_incr(1)*ddir
             do is = 2, msc -1 
               ETOT_SQ_WK = ETOT_SQ_WK + ONEHALF*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
             ETOT_SQ_WK = ETOT_SQ_WK + ONEHALF*tmp(msc) * ds_incr(msc)*ddir
           end do
           !tmp = SQRT(WK(:,IP))
           !ETOT_SQ_WK  = DINTSPEC_Y(IP,ACLOC,tmp) 
!
! tail factors ...
!
           DS          = SPSIG(MSC) - SPSIG(MSC-1)

           ATAIL       = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,1) * DDIR * DS
           ETAIL       = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,2) * DDIR * DS
           ESIGTAIL    = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,3) * DDIR * DS
!
! tail factors ...
!
           ACTOT       = ACTOT        + PTAIL(5)  * ATAIL 
           ETOT        = ETOT         + PTAIL(6)  * ETAIL 
           ETOT_SPSIG  = ETOT_SPSIG   + PTAIL(7)  * ETAIL 
           ETOT_ISQ_WK = ETOT_ISQ_WK  + PTAIL(5)  * ETAIL / (SQRT(WK(MSC,IP)))
           ETOT_SQ_WK  = ETOT_SQ_WK   + PTAIL(5)  * ETAIL * (SQRT(WK(MSC,IP)))
           ETOT_WK     = ETOT_WK      + PTAIL(8)  * ETAIL * WK(MSC,IP)
!
! integral parameters ...
!
           HS          = MAX(ZERO,4.*SQRT(ETOT))
           IF (LMONO_OUT) HS = HS / SQRT(2.)
           SME01       = ETOT_SPSIG / ETOT
           KME01       = ETOT_WK / ETOT
           SME10       = ETOT / ACTOT
           KMWAM       = (ETOT/ETOT_ISQ_WK)**2
           KMWAM2      = (ETOT_SQ_WK/ETOT)**2

         else
!
! no or too less energy ...
!
           ETOT        = ZERO
           HS          = ZERO 
           SME01       = ZERO 
           KME01       = 10.0_rkind
           SME10       = ZERO 
           KMWAM       = 10.0_rkind
           KMWAM2      = 10.0_rkind

         end if 
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_PARAMETER_BDCONS(ACLOC,HS,TM01,TM02)

      USE DATAPOOL
      IMPLICIT NONE

      REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)
      REAL(rkind), INTENT(OUT)   :: HS,TM01,TM02

      INTEGER             :: ID, IS

      REAL(rkind)                :: Y(MSC)
      REAL(rkind)                :: DS, ETAIL
      REAL(rkind)                :: OMEG2, EAD, ETOT
      REAL(rkind)                :: EFTAIL,PPTAIL,EFTOT,EPTAIL
      REAL(rkind)                :: EHFR,AHFR,APTAIL,EPTOT,APTOT
      REAL(rkind)                :: tmp(msc),actmp(msc)
!
! total energy ...
!
      Y = ZERO
      tmp = ZERO
      actmp = ZERO
      ETOT = ZERO
      do id = 1, mdc
        tmp(:) = acloc(:,id) * spsig
        ETOT = ETOT + tmp(1) * ONEHALF * ds_incr(1)*ddir
        do is = 2, msc
          ETOT = ETOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
        end do
        ETOT = ETOT + ONEHALF * tmp(msc) * ds_incr(msc)*ddir
      end do

      IF (ETOT .GT. THR) THEN
!
! tail ratios
!
         DS    = SPSIG(MSC) - SPSIG(MSC-1)
         ETAIL = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,2) * DDIR * DS
         ETOT  = ETOT + PTAIL(6) * ETAIL

         HS = 4*SQRT(ETOT)

         APTOT = ZERO
         EPTOT = ZERO
         PPTAIL = PTAIL(1)
         APTAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))
         PPTAIL = PTAIL(1) - ONE
         EPTAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))

         DO ID = 1, MDC
           DO IS = 1, MSC
             APTOT = APTOT + SPSIG(IS)    * ACLOC(IS,ID)
             EPTOT = EPTOT + SIGPOW(IS,2) * ACLOC(IS,ID)
           ENDDO
         ENDDO

         APTOT = APTOT * FRINTF
         EPTOT = EPTOT * FRINTF

         IF (MSC .GT. 3  .AND. .NOT. LSIGMAX) THEN
           DO ID = 1, MDC
           AHFR  = SPSIG(MSC) * ACLOC(MSC,ID)
           APTOT = APTOT + APTAIL * AHFR
           EHFR  = SPSIG(MSC) * AHFR
           EPTOT = EPTOT + EPTAIL * EHFR
           ENDDO
         ENDIF

         IF (EPTOT .GT. ZERO) THEN
            TM01 = PI2 * APTOT / EPTOT
         ELSE
            TM01 = ZERO
         END IF

         ETOT  = ZERO
         EFTOT = ZERO
         PPTAIL = PTAIL(1) - ONE
         ETAIL  = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))
         PPTAIL = PTAIL(1) - 3.
         EFTAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))
         DO ID=1, MDC
            DO IS = 1, MSC
              EAD  = SIGPOW(IS,2) * ACLOC(IS,ID) * FRINTF
              OMEG2 = SIGPOW(IS,2)
              ETOT  = ETOT + EAD
              EFTOT = EFTOT + EAD * OMEG2
            ENDDO
            IF (MSC .GT. 3  .AND. .NOT. LSIGMAX) THEN
              EAD  = SIGPOW(MSC,2) * ACLOC(MSC,ID)
              ETOT  = ETOT  + ETAIL * EAD
              EFTOT = EFTOT + EFTAIL * OMEG2 * EAD
            ENDIF
         ENDDO
         IF (EFTOT .GT. sqrt(verysmall)) THEN
           TM02 = PI2 * SQRT(ETOT/EFTOT)
         ELSE
           TM02 = ZERO
         END IF

      ELSE
        HS = ZERO
        TM01 = ZERO
        TM02 = ZERO
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_PARAMETER(IP,ACLOC,ISMAX,HS,TM01,TM02,TM10,KLM,WLM)

      USE DATAPOOL
      IMPLICIT NONE
!2do ... rewrite this integration ...
      INTEGER, INTENT(IN) :: IP,ISMAX

      REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)
      REAL(rkind), INTENT(OUT)   :: HS,TM01,TM02,KLM,WLM,TM10

      INTEGER             :: ID, IS

      REAL(rkind)                :: Y(MSC)
      REAL(rkind)                :: DS, ETAIL
      REAL(rkind)                :: OMEG2,OMEG,EAD,UXD, ETOT
      REAL(rkind)                :: EFTAIL,PPTAIL,EFTOT,EPTAIL
      REAL(rkind)                :: EHFR,AHFR,APTAIL,EPTOT,APTOT
      REAL(rkind)                :: SKK, CKTAIL, ETOT1, SIG22, EKTOT, CETAIL
      REAL(rkind)                :: tmp(msc),actmp(msc)
!
! total energy ...
!
      Y = ZERO
      tmp = ZERO
      actmp = ZERO
      ETOT = ZERO
      do id = 1, mdc
        tmp(:) = acloc(:,id) * spsig
        ETOT = ETOT + tmp(1) * ONEHALF * ds_incr(1)*ddir
        do is = 2, msc
          ETOT = ETOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
        end do
        ETOT = ETOT + ONEHALF * tmp(msc) * ds_incr(msc)*ddir
      end do

      IF (ETOT .GT. THR) THEN
!
! tail ratios
!
         DS    = SPSIG(MSC) - SPSIG(MSC-1)
         ETAIL = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,2) * DDIR * DS
         ETOT  = ETOT + PTAIL(6) * ETAIL

         HS = 4*SQRT(ETOT)

         APTOT = ZERO
         EPTOT = ZERO
         PPTAIL = PTAIL(1)
         APTAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))
         PPTAIL = PTAIL(1) - ONE
         EPTAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))

         DO ID = 1, MDC
           DO IS = 1, ISMAX
             APTOT = APTOT + SPSIG(IS)    * ACLOC(IS,ID)
             EPTOT = EPTOT + SIGPOW(IS,2) * ACLOC(IS,ID)
           ENDDO
         ENDDO

         APTOT = APTOT * FRINTF
         EPTOT = EPTOT * FRINTF

         IF (MSC .GT. 3  .AND. .NOT. LSIGMAX) THEN
           DO ID = 1, MDC
           AHFR  = SPSIG(MSC) * ACLOC(MSC,ID)
           APTOT = APTOT + APTAIL * AHFR
           EHFR  = SPSIG(MSC) * AHFR
           EPTOT = EPTOT + EPTAIL * EHFR
           ENDDO
         ENDIF

         IF (EPTOT .GT. ZERO) THEN
            TM01 = PI2 * APTOT / EPTOT
         ELSE
            TM01 = ZERO
         END IF

         ETOT  = ZERO
         EFTOT = ZERO
         PPTAIL = PTAIL(1) - ONE
         ETAIL  = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))
         PPTAIL = PTAIL(1) - 3.
         EFTAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))
         DO ID=1, MDC
            IF (LSECU .OR. LSTCU) THEN
              UXD  = CURTXY(IP,1)*COSTH(ID) + CURTXY(IP,2)*SINTH(ID)
            ENDIF
            DO IS = 1, ISMAX
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
            IF (MSC .GT. 3  .AND. .NOT. LSIGMAX) THEN
              EAD  = SIGPOW(MSC,2) * ACLOC(MSC,ID)
              ETOT  = ETOT  + ETAIL * EAD
              EFTOT = EFTOT + EFTAIL * OMEG2 * EAD
            ENDIF
         ENDDO
         IF (EFTOT .GT. sqrt(verysmall)) THEN
           TM02 = PI2 * SQRT(ETOT/EFTOT)
         ELSE
           TM02 = ZERO
         END IF

         ETOT1 = ZERO
         EKTOT = ZERO
!
! tail ratios same
!
         PPTAIL = PTAIL(1) - ONE
         CETAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))
         PPTAIL = PTAIL(1) - ONE - 2*ONE
         CKTAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))

         DO IS = 1, ISMAX
           SIG22 = SIGPOW(IS,2)
           SKK  = SIG22 * WK(IS,IP)
           DO ID = 1, MDC
             ETOT1 = ETOT1 + SIG22 * ACLOC(IS,ID)
             EKTOT = EKTOT + SKK * ACLOC(IS,ID)
           ENDDO
         ENDDO
         ETOT1 = FRINTF * ETOT1
         EKTOT = FRINTF * EKTOT

         IF (MSC .GT. 3) THEN
            DO ID=1,MDC
              ETOT1 = ETOT1 + CETAIL * SIG22 * ACLOC(MSC,ID)
              EKTOT = EKTOT + CKTAIL * SKK * ACLOC(MSC,ID)
            ENDDO
         ENDIF

         IF (ETOT1.GT.VERYSMALL.AND.EKTOT.GT.VERYSMALL) THEN
            WLM = PI2 * (ETOT1/EKTOT)
            KLM = PI2/WLM
         ELSE
            KLM = 10.0_rkind
            WLM = ZERO
         ENDIF

         APTOT = 0.
         EPTOT = 0.
         DO ID=1, MDC
           DO IS=1,ISMAX
             APTOT = APTOT + SPSIG(IS) * ACLOC(IS,ID)
             EPTOT = EPTOT + SIGPOW(IS,2) * ACLOC(IS,ID)
           ENDDO
         ENDDO
         APTOT = APTOT * FRINTF
         EPTOT = EPTOT * FRINTF
         IF (MSC .GT. 3) THEN
           PPTAIL = PTAIL(1)
           APTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
           PPTAIL = PTAIL(1) - 1.
           EPTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
           DO ID = 1, MDC
             AHFR = SPSIG(MSC) * ACLOC(MSC,ID)
             APTOT = APTOT + APTAIL * AHFR
             EHFR = SPSIG(MSC) * AHFR
             EPTOT = EPTOT + EPTAIL * EHFR
           ENDDO
         ENDIF
         TM10 = 2.*PI * APTOT / EPTOT

      ELSE

           HS = ZERO
           TM01 = ZERO
           TM02 = ZERO
           TM10 = ZERO
           KLM  = 10.0_rkind
           WLM  = ZERO

      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_WAVE_PARAMETER_LOC(ACLOC,CURTXYLOC,DEPLOC,WKLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)

         USE DATAPOOL
         IMPLICIT NONE
 
         REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(IN)    :: WKLOC(MSC), DEPLOC
         REAL(rkind), INTENT(IN)    :: CURTXYLOC(2)
         REAL(rkind), INTENT(OUT)   :: SME01, SME10, KME01
         REAL(rkind), INTENT(OUT)   :: KMWAM, KMWAM2
         REAL(rkind), INTENT(OUT)   :: HS

         INTEGER             :: ID, IS

         REAL(rkind)                :: ACTOT, ETOT
         REAL(rkind)                :: ETOT_SPSIG
         REAL(rkind)                :: ETOT_WK
         REAL(rkind)                :: ETOT_ISQ_WK
         REAL(rkind)                :: ETOT_SQ_WK

         REAL(rkind)                :: Y(MSC)
         REAL(rkind)                :: DS, ATAIL, ETAIL, ESIGTAIL
         REAL(rkind)                :: tmp(msc)
!
! total energy ...
!
         ETOT = ZERO
         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig
           ETOT = ETOT + tmp(1) * ONEHALF * ds_incr(1)*ddir
           do is = 2, msc
             ETOT = ETOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
           end do
           ETOT = ETOT + ONEHALF * tmp(msc) * ds_incr(msc)*ddir
         end do
!
! if etot too small skip ...
!
         if (etot .gt. verysmall) then

           ACTOT = ZERO
           y = ONE/SPSIG
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ACTOT  = ACTOT + tmp(1) * ONEHALF * ds_incr(1)*ddir
             do is = 2, msc
               ACTOT = ACTOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ACTOT  = ACTOT + ONEHALF * tmp(msc) * ds_incr(msc)*ddir
           end do

           ETOT_SPSIG = ZERO
           y = SIGPOW(:,1)
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ETOT_SPSIG = ETOT_SPSIG + tmp(1) * ONEHALF * ds_incr(1)*ddir
             do is = 2, msc
               ETOT_SPSIG = ETOT_SPSIG + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ETOT_SPSIG = ETOT_SPSIG + ONEHALF * tmp(msc) * ds_incr(msc)*ddir
           end do

           ETOT_WK = ZERO
           y = WKLOC(:)
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ETOT_WK = ETOT_WK + ONEHALF * tmp(1) * ds_incr(1)*ddir
             do is = 2, msc
               ETOT_WK = ETOT_WK + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ETOT_WK = ETOT_WK + ONEHALF * tmp(msc) * ds_incr(msc)*ddir
           end do

           ETOT_ISQ_WK = ZERO
           y = ONE/SQRT(WKLOC)
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ETOT_ISQ_WK = ETOT_ISQ_WK + ONEHALF * tmp(1) * ds_incr(1)*ddir
             do is = 2, msc -1 
               ETOT_ISQ_WK = ETOT_ISQ_WK+ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ETOT_ISQ_WK = ETOT_ISQ_WK+ONEHALF*tmp(msc) * ds_incr(msc)*ddir
           end do

           ETOT_SQ_WK = ZERO
           y = SQRT(WKLOC)
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ETOT_SQ_WK = ETOT_SQ_WK + ONEHALF * tmp(1) * ds_incr(1)*ddir
             do is = 2, msc -1 
               ETOT_SQ_WK = ETOT_SQ_WK + ONEHALF*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
             ETOT_SQ_WK = ETOT_SQ_WK + ONEHALF*tmp(msc) * ds_incr(msc)*ddir
           end do
!
           DS          = SPSIG(MSC) - SPSIG(MSC-1)

           ATAIL       = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,1) * DDIR * DS
           ETAIL       = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,2) * DDIR * DS
           ESIGTAIL    = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,3) * DDIR * DS
!
! tail factors ...
!
           ACTOT       = ACTOT        + PTAIL(5)  * ATAIL 
           ETOT        = ETOT         + PTAIL(6)  * ETAIL 
           ETOT_SPSIG  = ETOT_SPSIG   + PTAIL(7)  * ETAIL 
           ETOT_ISQ_WK = ETOT_ISQ_WK  + PTAIL(5)  * ETAIL / SQRT(WKLOC(MSC))
           ETOT_SQ_WK  = ETOT_SQ_WK   + PTAIL(5)  * ETAIL * SQRT(WKLOC(MSC))
           ETOT_WK     = ETOT_WK      + PTAIL(8)  * ETAIL * WKLOC(MSC) 
!
! integral parameters ...
!
           HS          = MAX(ZERO,4.*SQRT(ETOT))
           SME01       = ETOT_SPSIG / ETOT
           KME01       = ETOT_WK / ETOT
           SME10       = ETOT / ACTOT
           KMWAM       = (ETOT/ETOT_ISQ_WK)**2
           KMWAM2      = (ETOT_SQ_WK/ETOT)**2

         else
!
! no or too less energy ...
!
           HS          = ZERO 
           SME01       = ZERO 
           KME01       = 10.0_rkind
           SME10       = ZERO 
           KMWAM       = 10.0_rkind
           KMWAM2      = 10.0_rkind

         end if 
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT,CALLFROM)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP
         REAL(rkind),    INTENT(IN) :: ACLOC(MSC,MDC)
         CHARACTER(len=*), INTENT(IN) :: CALLFROM

         REAL(rkind), INTENT(OUT)   :: UBOT, ORBITAL, BOTEXPER, TMBOT

         INTEGER             :: ID, IS

         REAL(rkind)                :: ETOT_SKD
         REAL(rkind)                :: ETOT_SKDSIG, TMP(MSC), Y(MSC)

!
! integrals ...
!
         IF (DEP(IP) .LT. DMIN) RETURN
         
         ETOT_SKD    = ZERO
         ETOT_SKDSIG = ZERO

         y = ONE/SINH(MIN(KDMAX,WK(:,IP)*DEP(IP)))**2

         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig * y
           ETOT_SKD  = ETOT_SKD + tmp(1) * ONEHALF * ds_incr(1)*ddir
           do is = 2, msc -1 
             ETOT_SKD = ETOT_SKD + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
           end do
           ETOT_SKD = ETOT_SKD + tmp(msc) * ONEHALF * ds_incr(msc)*ddir
         end do
 
         y =  SIGPOW(:,2)*ONE/SINH(MIN(KDMAX,WK(:,IP)*DEP(IP)))**2

         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig * y
           ETOT_SKDSIG = ETOT_SKDSIG + tmp(1) * ONEHALF * ds_incr(1)*ddir
           do is = 2, msc -1 
             ETOT_SKDSIG = ETOT_SKDSIG + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
           end do 
           ETOT_SKDSIG = ETOT_SKDSIG + tmp(msc) * ONEHALF * ds_incr(msc)*ddir
         end do
 
         IF (ETOT_SKD .gt. verysmall) THEN 
!
! integral parameters ...
!
           UBOT        = SQRT(ETOT_SKD)
           ORBITAL     = SQRT(2*ETOT_SKD)
           BOTEXPER    = SQRT(2*ETOT_SKDSIG)
           TMBOT       = PI2*SQRT(ETOT_SKDSIG/ETOT_SKD)
           
         ELSE 
!
! no or too less energy ...
! 
           UBOT        = ZERO 
           ORBITAL     = ZERO 
           BOTEXPER    = ZERO 
           TMBOT       = ZERO 

         ENDIF 
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_CURRENT_PARAMETER_LOC(ACLOC,CURTXYLOC,DEPLOC,WKLOC,UBOT,ORBITAL,BOTEXPER,TMBOT)
         USE DATAPOOL
         IMPLICIT NONE

         REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(IN)    :: WKLOC(MSC), DEPLOC
         REAL(rkind), INTENT(IN)    :: CURTXYLOC(2)

         REAL(rkind), INTENT(OUT)   :: UBOT, ORBITAL, BOTEXPER, TMBOT

         INTEGER             :: ID, IS

         LOGICAL             :: ISINF

         REAL(rkind)                :: ETOT_SKD
         REAL(rkind)                :: ETOT_SKDSIG, TMP(MSC), Y(MSC)
!
! integrals ...
!
         IF (DEPLOC .LT. DMIN) RETURN

         ETOT_SKD    = ZERO
         ETOT_SKDSIG = ZERO

         y = ONE/SINH(MIN(KDMAX,WKLOC*DEPLOC))**2
         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig * y
           ETOT_SKD  = ETOT_SKD + tmp(1) * ONEHALF * ds_incr(1)*ddir
           do is = 2, msc -1
             ETOT_SKD = ETOT_SKD + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
           end do
           ETOT_SKD = ETOT_SKD + tmp(msc) * ONEHALF * ds_incr(msc)*ddir
         end do

         y =  SIGPOW(:,2)*ONE/SINH(MIN(KDMAX,WKLOC*DEPLOC))**2
         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig * y
           ETOT_SKDSIG = ETOT_SKDSIG + tmp(1) * ONEHALF * ds_incr(1)*ddir
           do is = 2, msc -1
             ETOT_SKDSIG = ETOT_SKDSIG + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
           end do
           ETOT_SKDSIG = tmp(msc) * ONEHALF * ds_incr(msc)*ddir
         end do

         IF (ETOT_SKD .gt. verysmall) THEN 
!
! integral parameters ...
!
           UBOT        = SQRT(ETOT_SKD)
           ORBITAL     = SQRT(2*ETOT_SKD)
           BOTEXPER    = SQRT(2*ETOT_SKDSIG)
           TMBOT       = PI2*SQRT(ETOT_SKDSIG/ETOT_SKD)
           
         ELSE 
!
! no or too less energy ...
! 
           UBOT        = ZERO 
           ORBITAL     = ZERO 
           BOTEXPER    = ZERO 
           TMBOT       = ZERO 

         ENDIF 

         !WRITE(*,'(9F15.4)') DEPLOC,SUM(ACLOC),CURTXYLOC,SUM(WKLOC), ETOT_SKD, ETOT_SKDSIG 
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE URSELL_NUMBER(HS,SME,DEPTH,URSELL)
         USE DATAPOOL, ONLY : G9, DMIN, verysmall, rkind, ONE, TWO, ZERO

         IMPLICIT NONE

         REAL(rkind), INTENT(IN)    :: HS, SME, DEPTH
         REAL(rkind), INTENT(OUT)   :: URSELL 

         IF (DEPTH .GT. DMIN .AND. SME .GT. verysmall .AND. HS .GT. verysmall) THEN
           URSELL = (G9 * HS)/(TWO*SQRT(TWO)*SME**2*DEPTH**2)
         ELSE
           URSELL = ZERO
         END IF 
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WINDSEASWELLSEP( IP, ACLOC, TM_W, CGP_W, CP_W, TP_W, LP_W, HS_W, KP_W )
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)           :: IP
         REAL(rkind)   , INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind)   , INTENT(OUT)   :: TM_W, CGP_W, CP_W, TP_W, LP_W, HS_W, KP_W

         INTEGER                :: ID, IS

         REAL(rkind)                   :: ETOT
         REAL(rkind)                   :: VEC2RAD, EFTOT, OMEG, WINDTH
         REAL(rkind)                   :: EAD, DS, EHFR, EFTAIL, ETAIL, PPTAIL, ACWIND(MSC,MDC)
         REAL(rkind)                   :: UXD, ETOTF3, ETOTF4, FP, WN_W, WVC, WKDEP_W

         WINDTH = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))

         ETOT = ZERO
         EFTAIL = ONE / (PTAIL(1)-ONE)

         DO ID  = 1, MDC            ! Calculate wind sea energy ... weak criterion
           DO IS = 1, MSC
             WVC = MyREAL(SPSIG(IS)/WK(IS,IP))
             IF (  1.2*UFRIC(IP)*COS(SPDIR(ID)-WINDTH)*(28./WVC) .LT. ONE) THEN
               ACWIND(IS,ID) = ZERO   ! Swell
             ELSE
               ACWIND(IS,ID) = ACLOC(IS,ID)  ! Wind Sea
             END IF
           END DO
           DO IS = 2, MSC
              DS = SPSIG(IS) - SPSIG(IS-1)
              EAD = ONEHALF*(SPSIG(IS)*ACWIND(IS,ID)+SPSIG(IS-1)*ACWIND(IS-1,ID))*DS*DDIR
              ETOT = ETOT + EAD
           END DO
           IF (MSC > 3) THEN
             EHFR = ACLOC(MSC,ID) * SPSIG(MSC)
             ETOT = ETOT + DDIR * EHFR * SPSIG(MSC) * EFTAIL
           ENDIF
         END DO

         IF (ETOT .LT. VERYSMALL) RETURN

         IF (ETOT > ZERO) THEN
           HS_W = 4.0*SQRT(ETOT)
         ELSE
           HS_W = ZERO
         END IF

         ETOTF3 = ZERO
         ETOTF4 = ZERO
         DO IS = 1, MSC
           DO ID = 1, MDC
             ETOTF3 = ETOTF3 + SPSIG(IS) * ACLOC(IS,ID)**4 * DDIR * DS_BAND(IS)
             ETOTF4 = ETOTF4 +             ACLOC(IS,ID)**4 * DDIR * DS_BAND(IS)
           END DO
         END DO
         IF(ETOTF4 .GT. VERYSMALL) THEN
            FP = ETOTF3/ETOTF4
            TP_W = ONE/FP/PI2
            !CALL WAVEKCG(DEP(IP), FP, WN_W, CP_W, KP_W, CGP_W)
            CALL ALL_FROM_TABLE(FP,DEP(IP),KP_W,CGP_W,WKDEP_W,WN_W,CP_W)
            LP_W  = PI2/KP_w
         ELSE
            FP    = ZERO 
            CP_W  = ZERO 
            KP_W  = 10.0_rkind
            CGP_W = ZERO 
            LP_W  = ZERO 
         END IF

         ETOT = ZERO
         EFTOT = ZERO
         PPTAIL = PTAIL(1) - ONE
         ETAIL  = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))
         PPTAIL = PTAIL(1) - 2.
         EFTAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))
         DO ID = 1, MDC
            UXD = CURTXY(IP,1)*COSTH(ID) + CURTXY(IP,2)*SINTH(ID)
            DO IS = 1, MSC
              OMEG = SPSIG(IS) + WK(IS,IP) * UXD
              EAD = FRINTF * SIGPOW(IS,2) * ACWIND(IS,ID)
              ETOT = ETOT + EAD
              EFTOT = EFTOT + EAD * OMEG
            ENDDO
            IF (MSC .GT. 3) THEN
              EAD = SIGPOW(MSC,2) * ACWIND(MSC,ID)
              ETOT = ETOT + ETAIL * EAD
             EFTOT = EFTOT + EFTAIL * OMEG * EAD
            ENDIF
         ENDDO
         IF (EFTOT.GT.ZERO) THEN
            TM_W = PI2 * ETOT / EFTOT
         ELSE
            TM_W = ZERO
         ENDIF
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PEAK_PARAMETER(IP,ACLOC,ISMAX,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)           :: IP, ISMAX
         REAL(rkind), INTENT(IN)       :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(OUT)      :: KPP, FPP, CPP, WNPP, CGPP, TPP, LPP, PEAKDSPR,PEAKDM,DPEAK
         REAL(rkind), INTENT(OUT)      :: TPPD,KPPD,CGPD,CPPD

         INTEGER                       :: IS, ID, IDIRM, ISIGMP
         REAL(rkind)                   :: HQUOT, HQUOTP, ETOTF3, ETOTF4, ETOTC4, ETOTS4, PEAKFF,WKDEPD,WNPD
         REAL(rkind)                   :: DEG, VEC2DEG, MAXAC, E1, E2, DS, EAD, ETOTT, WKDEPP
!
! Peak period continues version... Taken from Thesis Henriques Alves ... correct citation is given there ... :)
!
       MAXAC = MAXVAL(ACLOC)

       IF (MAXAC .gt. VERYSMALL .AND.  DEP(IP) .GT. DMIN) THEN

         ETOTF3 = ZERO
         ETOTF4 = ZERO
         ETOTC4 = ZERO
         ETOTS4 = ZERO

         DO IS = 1, MSC
           DO ID = 1, MDC
             HQUOT  = ACLOC(IS,ID)/MAXAC
             HQUOTP = HQUOT**4
             ETOTF3 = ETOTF3 + SPSIG(IS) * HQUOTP * DS_BAND(IS)
             ETOTF4 = ETOTF4 +             HQUOTP * DS_BAND(IS)
             ETOTC4 = ETOTC4 + COSTH(ID) * HQUOTP * DS_BAND(IS)
             ETOTS4 = ETOTS4 + SINTH(ID) * HQUOTP * DS_BAND(IS)
           END DO
         END DO

         IF(ETOTF4 .GT. VERYSMALL .AND. ETOTF4 .GT. VERYSMALL) THEN

           FPP    = ETOTF3/ETOTF4
           !CALL WAVEKCG(DEP(IP), FPP, WNPP, CPP, KPP, CGPP)
           CALL ALL_FROM_TABLE(FPP,DEP(IP),KPP,CGPP,WKDEPP,WNPP,CPP)
           TPP    = PI2/FPP
           LPP    = PI2/KPP
           PEAKDM = VEC2DEG (ETOTC4, ETOTS4)
           CALL DEG2NAUT(PEAKDM,DEG,LNAUTOUT)
           PEAKDM = DEG
           PEAKFF = MIN(ONE,SQRT(ETOTC4*ETOTC4+ETOTS4*ETOTS4)/ETOTF4)
           PEAKDSPR = SQRT(MAX(ZERO,2.-2.*PEAKFF)) * 180./PI

         ELSE 

           FPP = ZERO
           KPP = 10.0_rkind
           CGPP = ZERO 
           WKDEPP = ZERO
           WNPP = ZERO
           CPP = ZERO
           TPP = ZERO
           LPP = ZERO
           PEAKDM = ZERO
           PEAKFF = ZERO
           PEAKDSPR = ZERO 

         END IF

         DPEAK = 1
         ETOTT = ZERO
         IDIRM = -1
         DO ID = 1, MDC
            EAD = ZERO
            DO IS = 2, ISMAX
               DS = SPSIG(IS)-SPSIG(IS-1)
               E1 = SPSIG(IS-1)*ACLOC(IS-1,ID)
               E2 = SPSIG(IS)*ACLOC(IS,ID)
               EAD = EAD + DS*(E1+E2)
            END DO
            IF (EAD .GT. ETOTT) THEN
               ETOTT = EAD
               IDIRM = ID
            END IF
         END DO
         IF (IDIRM .GT. 0) THEN
           DPEAK    = SPDIR(IDIRM)
         ELSE
           DPEAK    = ZERO
         END IF

       ELSE

         FPP = ZERO
         KPP = 10.0_rkind
         CGPP = ZERO
         WKDEPP = ZERO
         WNPP = ZERO
         CPP = ZERO
         TPP = ZERO
         LPP = ZERO
         PEAKDM = ZERO
         PEAKFF = ZERO
         PEAKDSPR = ZERO
         DPEAK = ZERO

       END IF


       ETOTT = ZERO
       ISIGMP = -1
       DO IS = 1, MSC
         EAD = ZERO
         DO ID = 1, MDC
            EAD = EAD + SPSIG(IS)*ACLOC(IS,ID)*DDIR
         ENDDO
         IF (EAD > ETOTT) THEN
           ETOTT = EAD
           ISIGMP = IS
          END IF
       END DO
       IF (ISIGMP > 0) THEN
          TPPD = ONE/(SPSIG(ISIGMP)/PI2)
          !CALL WAVEKCG(DEP(IP), SPSIG(ISIGMP), CPPD, KPPD, CGPD)
          CALL ALL_FROM_TABLE(SPSIG(ISIGMP),DEP(IP),KPPD,CGPD,WKDEPD,WNPD,CPPD)
       ELSE
          TPPD = ZERO
          CPPD  = ZERO
          KPPD  = ZERO
          CGPD  = ZERO
       END IF

       !WRITE(*,'(11F15.4)') FPP, KPP, CGPP, WKDEPP, WNPP, CPP, TPP, LPP, PEAKDM, PEAKFF, PEAKDSPR 
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PEAK_PARAMETER_LOC(ACLOC,DEPLOC,ISMAX,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: ISMAX
         REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(IN)    :: DEPLOC


         REAL(rkind)  , INTENT(OUT)    :: KPP,CPP,WNPP,CGPP,TPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD,FPP

         INTEGER                :: IS, ID, IDIRM, ISIGMP
         REAL(rkind)            :: HQUOT, HQUOTP, ETOTF3, ETOTF4, ETOTC4, ETOTS4, PEAKFF, WKDEPP
         REAL(rkind)            :: DEG, VEC2DEG, MAXAC, E1, E2, DS, EAD, ETOTT,WKDEPD,WNPD
!
! Peak period continues version... Taken from Thesis Henriques Alves ... correct citation is given there ... :)
!
       MAXAC = MAXVAL(ACLOC)

       IF (MAXAC .gt. VERYSMALL .AND. DEPLOC .GT. DMIN) THEN
         ETOTF3 = ZERO
         ETOTF4 = ZERO
         ETOTC4 = ZERO
         ETOTS4 = ZERO
         DO IS = 1, MSC
           DO ID = 1, MDC
             HQUOT  = ACLOC(IS,ID)/MAXAC
             HQUOTP = HQUOT**4
             ETOTF3 = ETOTF3 + SPSIG(IS) * HQUOTP * DS_BAND(IS)
             ETOTF4 = ETOTF4 +             HQUOTP * DS_BAND(IS)
             ETOTC4 = ETOTC4 + COSTH(ID) * HQUOTP * DS_BAND(IS)
             ETOTS4 = ETOTS4 + SINTH(ID) * HQUOTP * DS_BAND(IS)
           END DO
         END DO
         IF(ETOTF4 .GT. VERYSMALL .AND. ETOTF4 .GT. VERYSMALL) THEN
           FPP    = ETOTF3/ETOTF4
           !CALL WAVEKCG(DEPLOC, FPP, WNPP, CPP, KPP, CGPP)
           CALL ALL_FROM_TABLE(FPP,DEPLOC,KPP,CGPP,WKDEPP,WNPP,CPP) 
           PEAKDM = VEC2DEG (ETOTC4, ETOTS4)
           TPP    = PI2/FPP
           LPP    = PI2/KPP
           CALL DEG2NAUT(PEAKDM,DEG,LNAUTOUT)
           PEAKDM = DEG
           PEAKFF = MIN(ONE,SQRT(MAX(ZERO,ETOTC4*ETOTC4+ETOTS4*ETOTS4))/ETOTF4)
           PEAKDSPR = SQRT(MAX(ZERO,2.-2.*PEAKFF)) * 180./PI
         ELSE
           FPP = ZERO
           KPP = 10.0_rkind
           CGPP = ZERO
           WKDEPP = ZERO
           WNPP = ZERO
           CPP = ZERO
           TPP = ZERO
           LPP = ZERO
           PEAKDM = ZERO
           PEAKFF = ZERO
           PEAKDSPR = ZERO
         END IF
         DPEAK = 1
         ETOTT = ZERO
         IDIRM = -1
         DO ID = 1, MDC
            EAD = ZERO
            DO IS = 2, ISMAX
               DS = SPSIG(IS)-SPSIG(IS-1)
               E1 = SPSIG(IS-1)*ACLOC(IS-1,ID)
               E2 = SPSIG(IS)*ACLOC(IS,ID)
               EAD = EAD + DS*(E1+E2)
            END DO
            IF (EAD .GT. ETOTT) THEN
               ETOTT = EAD
               IDIRM = ID
            END IF
         END DO
         IF (IDIRM .GT. 0) THEN
           DPEAK    = SPDIR(IDIRM) * RADDEG
           CALL DEG2NAUT(DPEAK,DEG,LNAUTOUT)
           DPEAK = DEG
         ELSE
           DPEAK = ZERO
         END IF
       ELSE
         FPP = ZERO
         KPP = 10.0_rkind
         CGPP = ZERO
         WKDEPP = ZERO
         WNPP = ZERO
         CPP = ZERO
         TPP = ZERO
         LPP = ZERO
         PEAKDM = ZERO
         PEAKFF = ZERO
         PEAKDSPR = ZERO
         DPEAK = ZERO
       END IF

       ETOTT = ZERO
       ISIGMP = -1
       DO IS = 1, MSC
         EAD = ZERO
         DO ID = 1, MDC
            EAD = EAD + SPSIG(IS)*ACLOC(IS,ID)*DDIR
         ENDDO
         IF (EAD > ETOTT) THEN
           ETOTT = EAD
           ISIGMP = IS
         END IF
       END DO
       IF (ISIGMP > 0) THEN
          TPPD = ONE/(SPSIG(ISIGMP)/PI2)
          !CALL WAVEKCG(DEPLOC, SPSIG(ISIGMP), CPPD, KPPD, CGPD)
          CALL ALL_FROM_TABLE(SPSIG(ISIGMP),DEPLOC,KPPD,CGPD,WKDEPD,WNPD,CPPD)
       ELSE
          TPPD = ZERO
          CPPD  = ZERO
          KPPD  = ZERO
          CGPD  = ZERO
       END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,ISMAX,ETOTS,ETOTC,DM,DSPR)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP, ISMAX
         REAL(rkind),    INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind),    INTENT(OUT)   :: ETOTS, ETOTC, DM, DSPR

         INTEGER                       :: IS, ID
         REAL(rkind)                   :: DS, EDI, EAD, ETOT1, EHFR
         REAL(rkind)                   :: EFTAIL, VEC2DEG, DEG, FF

         ETOTC = ZERO
         ETOTS = ZERO
         ETOT1 = ZERO

         EFTAIL = ONE / (PTAIL(1)-ONE)

         DO ID = 1, MDC
           EAD = ZERO
             DO  IS = 2, ISMAX 
               DS  = SPSIG(IS)-SPSIG(IS-1)
               EDI = ONEHALF*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS
               EAD = EAD + EDI
            END DO
            IF (MSC .GT. 3) THEN
              EHFR = ACLOC(MSC,ID) * SPSIG(MSC)
              EAD = EAD + EHFR * SPSIG(MSC) * EFTAIL
            ENDIF
            EAD = EAD * DDIR
            ETOT1 = ETOT1 + EAD
            ETOTC  = ETOTC + EAD * COSTH(ID)
            ETOTS  = ETOTS + EAD * SINTH(ID)
         END DO

         IF (ETOT1 .GT. verysmall) THEN
           DM    = VEC2DEG (ETOTC, ETOTS)
           CALL DEG2NAUT(DM,DEG,LNAUTOUT)
           DM = DEG
           FF = MIN (ONE, SQRT(MAX(ZERO,ETOTC*ETOTC+ETOTS*ETOTS))/ETOT1)
           DSPR = SQRT(MAX(ZERO,2.-2.*FF)) * 180./PI
         ELSE
           FF   = ZERO
           DM   = ZERO
           DSPR = ZERO
         END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_DIRECTION_AND_SPREAD_LOC(ACLOC,ISMAX,ETOTS,ETOTC,DM,DSPR)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: ISMAX
         REAL(rkind),    INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind),    INTENT(OUT)   :: ETOTS, ETOTC, DM, DSPR

         INTEGER                :: IS, ID

         REAL(rkind)                   :: DS, EDI, EAD, ETOT1, EHFR
         REAL(rkind)                  :: EFTAIL, VEC2DEG, DEG, FF

         ETOTC = ZERO
         ETOTS = ZERO
         ETOT1  = ZERO

         EFTAIL = ONE / (PTAIL(1)-ONE)

         DO ID = 1, MDC
           EAD = ZERO
             DO  IS = 2, ISMAX 
               DS  = SPSIG(IS)-SPSIG(IS-1)
               EDI = ONEHALF*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS
               EAD = EAD + EDI
            END DO
            IF (MSC .GT. 3) THEN
              EHFR = ACLOC(MSC,ID) * SPSIG(MSC)
              EAD = EAD + EHFR * SPSIG(MSC) * EFTAIL
            ENDIF
            EAD = EAD * DDIR
            ETOT1 = ETOT1 + EAD
            ETOTC  = ETOTC + EAD * COSTH(ID)
            ETOTS  = ETOTS + EAD * SINTH(ID)
         END DO

         IF (ETOT1 .GT. verysmall ) THEN
           DM    = VEC2DEG (ETOTC, ETOTS)
           CALL DEG2NAUT(DM,DEG,LNAUTOUT)
           DM = DEG
           FF = MIN (ONE, SQRT(MAX(ZERO,ETOTC*ETOTC+ETOTS*ETOTS))/ETOT1)
           DSPR = SQRT(MAX(ZERO,2.-2.*FF)) * 180./PI
         ELSE
           FF   = ZERO
           DM   = ZERO
           DSPR = ZERO
         END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_PARAMETER_LOC(ACLOC,CURTXYLOC,DEPLOC,WKLOC,ISMAX,HS,TM01,TM02,TM10,KLM,WLM)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ISMAX
      REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(IN)    :: WKLOC(MSC), DEPLOC
         REAL(rkind), INTENT(IN)    :: CURTXYLOC(2)


         REAL(rkind), INTENT(OUT)   :: HS,TM01,TM02,KLM,WLM,TM10

         INTEGER             :: ID, IS

         REAL(rkind)                :: DS,ETAIL
         REAL(rkind)                :: OMEG2,OMEG,EAD,UXD, ETOT
         REAL(rkind)                :: EFTAIL,PPTAIL,EFTOT,EPTAIL
         REAL(rkind)                :: EHFR,AHFR,APTAIL,EPTOT,APTOT
         REAL(rkind)                :: SKK, CKTAIL, ETOT1, SIG22, EKTOT, CETAIL
         REAL(rkind)                :: tmp(msc)
!
! total energy ...
!
         ETOT = ZERO
         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig
           ETOT = ETOT + tmp(1) * ONEHALF * ds_incr(1)*ddir
           do is = 2, msc
             ETOT = ETOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
           end do
           ETOT = ETOT + ONEHALF * tmp(msc) * ds_incr(msc)*ddir
         end do

         IF (ETOT .GT. verysmall) THEN
!
! tail ratios same as in swan ...
!
         DS    = SPSIG(MSC) - SPSIG(MSC-1)
         ETAIL = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,2) * DDIR * DS
         ETOT  = ETOT + PTAIL(6) * ETAIL

         HS = 4*SQRT(ETOT)

         APTOT = ZERO
         EPTOT = ZERO
!
! tail ratios same as in swan ...
!
         PPTAIL = PTAIL(1)
         APTAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))
         PPTAIL = PTAIL(1) - ONE
         EPTAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))

         DO ID = 1, MDC
           DO IS = 1, ISMAX
             APTOT = APTOT + SPSIG(IS)    * ACLOC(IS,ID)
             EPTOT = EPTOT + SIGPOW(IS,2) * ACLOC(IS,ID)
           ENDDO
         ENDDO

         APTOT = APTOT * FRINTF
         EPTOT = EPTOT * FRINTF

         IF (MSC .GT. 3  .AND. .NOT. LSIGMAX) THEN
           DO ID = 1, MDC
             AHFR  = SPSIG(MSC) * ACLOC(MSC,ID)
             APTOT = APTOT + APTAIL * AHFR
             EHFR  = SPSIG(MSC) * AHFR
             EPTOT = EPTOT + EPTAIL * EHFR
           ENDDO
         ENDIF

         IF (EPTOT .GT. VERYSMALL) THEN
            TM01 = PI2 * APTOT / EPTOT
         ELSE
            TM01 = ZERO
         END IF

         ETOT  = ZERO
         EFTOT = ZERO

         PPTAIL = PTAIL(1) - ONE
         ETAIL  = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))
         PPTAIL = PTAIL(1) - 3.
         EFTAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))
!
! tail ratios same as in swan ...
!
         DO ID=1, MDC
            IF (LSECU .OR. LSTCU) THEN
              UXD  = CURTXYLOC(1)*COSTH(ID) + CURTXYLOC(2)*SINTH(ID)
            ENDIF
            DO IS = 1, ISMAX
              EAD  = SIGPOW(IS,2) * ACLOC(IS,ID) * FRINTF
              IF (LSECU .OR. LSTCU) THEN
                OMEG  = SPSIG(IS) + WKLOC(IS) * UXD
                OMEG2 = OMEG**2
              ELSE
               OMEG2 = SIGPOW(IS,2)
              ENDIF
              ETOT  = ETOT + EAD
              EFTOT = EFTOT + EAD * OMEG2
            ENDDO
            IF (MSC .GT. 3  .AND. .NOT. LSIGMAX) THEN
              EAD  = SIGPOW(MSC,2) * ACLOC(MSC,ID)
              ETOT  = ETOT  + ETAIL * EAD
              EFTOT = EFTOT + EFTAIL * OMEG2 * EAD
            ENDIF
         ENDDO
         IF (EFTOT .GT. verysmall .AND. ETOT .GT. verysmall) THEN
           TM02 = PI2 * SQRT(ETOT/EFTOT)
         ELSE
           TM02 = ZERO
         END IF

         ETOT1 = ZERO
         EKTOT = ZERO
!
! tail ratios
!
         PPTAIL = PTAIL(1) - ONE
         CETAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))
         PPTAIL = PTAIL(1) - ONE - 2.*ONE
         CKTAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))

         DO IS = 1, ISMAX
           SIG22 = SIGPOW(IS,2)
           SKK  = SIG22 * WKLOC(IS)
           DO ID = 1, MDC
             ETOT1 = ETOT1 + SIG22 * ACLOC(IS,ID)
             EKTOT = EKTOT + SKK * ACLOC(IS,ID)
           ENDDO
         ENDDO

         ETOT1 = FRINTF * ETOT1
         EKTOT = FRINTF * EKTOT

         IF (MSC .GT. 3) THEN
            DO ID=1,MDC
              ETOT1 = ETOT1 + CETAIL * SIG22 * ACLOC(MSC,ID)
              EKTOT = EKTOT + CKTAIL * SKK * ACLOC(MSC,ID)
            ENDDO
         ENDIF

         IF (ETOT1.GT.VERYSMALL.AND.EKTOT.GT.VERYSMALL) THEN
            WLM = PI2 * (ETOT1/EKTOT)
            KLM  = PI2/WLM
         ELSE
            KLM  = 10.0_rkind
            WLM = ZERO
         ENDIF

         APTOT = 0.
         EPTOT = 0.
         DO ID=1, MDC
           DO IS=1,ISMAX
             APTOT = APTOT + SPSIG(IS) * ACLOC(IS,ID)
             EPTOT = EPTOT + SIGPOW(IS,2) * ACLOC(IS,ID)
           ENDDO
         ENDDO
         APTOT = APTOT * FRINTF
         EPTOT = EPTOT * FRINTF
         IF (MSC .GT. 3) THEN
           PPTAIL = PTAIL(1)
           APTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
           PPTAIL = PTAIL(1) - 1.
           EPTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
           DO ID = 1, MDC
             AHFR = SPSIG(MSC) * ACLOC(MSC,ID)
             APTOT = APTOT + APTAIL * AHFR
             EHFR = SPSIG(MSC) * AHFR
             EPTOT = EPTOT + EPTAIL * EHFR
           ENDDO
         ENDIF
         TM10 = PI2 * APTOT / EPTOT

         ELSE

           HS   = ZERO
           TM01 = ZERO
           TM02 = ZERO
           TM10 = ZERO
           KLM  = 10.0_rkind
           WLM  = ZERO

         END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_FREQS(IP,ACLOC,SME01,SME10,ETOTWS,LWINDSEA)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP

         REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: SME01, SME10, ETOTWS

         LOGICAL, INTENT(IN) :: LWINDSEA(MSC,MDC)

         INTEGER             :: ID, IS

         REAL(rkind)                :: ACTOT, ETOT
         REAL(rkind)                :: ETOT_SPSIG

         REAL(rkind)                :: Y(MSC)
         REAL(rkind)                :: DS, ATAIL, ETAIL
         REAL(rkind)                :: tmp(msc)
!
! total energy ...
!
         ETOT = ZERO
         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig 
           ETOT = ETOT + ONEHALF * tmp(1) * ds_incr(1)*ddir
           do is = 2, msc
             IF (LWINDSEA(IS,ID)) ETOT = ETOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir 
           end do
           ETOT = ETOT + ONEHALF * tmp(msc) * ds_incr(msc)*ddir
         end do
!
! if etot too small skip ...
!
         if (etot .gt. verysmall) then
!
! integrals ... inlined ... for speed ...
!
           ACTOT = ZERO
           y = ONE/SPSIG
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ACTOT = ACTOT + ONEHALF * tmp(1) * ds_incr(1)*ddir
             do is = 2, msc
               IF (LWINDSEA(IS,ID)) ACTOT = ACTOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir 
             end do
             ACTOT = ACTOT + ONEHALF * tmp(msc) * ds_incr(msc)*ddir
           end do

           ETOT_SPSIG = ZERO
           y = SIGPOW(:,1) 
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ETOT_SPSIG = ETOT_SPSIG + ONEHALF * tmp(1) * ds_incr(1)*ddir
             do is = 2, msc
               IF (LWINDSEA(IS,ID)) ETOT_SPSIG = ETOT_SPSIG + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir 
             end do
             ETOT_SPSIG = ETOT_SPSIG + ONEHALF * tmp(msc) * ds_incr(msc)*ddir
           end do
!
! tail factors ...
!
           DS          = SPSIG(MSC) - SPSIG(MSC-1)

           ATAIL       = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,1) * DDIR * DS
           ETAIL       = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,2) * DDIR * DS
!
! tail factors ...
!
           ACTOT       = ACTOT        + PTAIL(5)  * ATAIL 
           ETOT        = ETOT         + PTAIL(6)  * ETAIL 
           ETOT_SPSIG  = ETOT_SPSIG   + PTAIL(7)  * ETAIL 
!
! integral parameters ...
!
           SME01       = ETOT_SPSIG / ETOT
           SME10       = ETOT / ACTOT
           ETOTWS      = ETOT

         else
!
! no or too less energy ...
!
           SME01       = ZERO 
           SME10       = ZERO 
           ETOTWS      = ZERO

         end if 
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_WAVEN(IP,ACLOC,KME01,KMWAM,KMWAM2)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP

         REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: KME01
         REAL(rkind), INTENT(OUT)   :: KMWAM, KMWAM2
         INTEGER             :: ID, IS


         REAL(rkind)                :: ACTOT, ETOT
         REAL(rkind)                :: ETOT_SPSIG
         REAL(rkind)                :: ETOT_WK
         REAL(rkind)                :: ETOT_ISQ_WK
         REAL(rkind)                :: ETOT_SQ_WK

         REAL(rkind)                :: Y(MSC), tmp(msc)
         REAL(rkind)                :: DS, ATAIL, ETAIL, ESIGTAIL
!         REAL(rkind)                :: dintspec, dintspec_y
!
! total energy ...
!
         !ETOT = DINTSPEC(IP,ACLOC)
         ETOT = ZERO
         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig
           do is = 2, msc
             ETOT = ETOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
           end do
         end do
!
! if etot too small skip ...
!
         if (etot .gt. verysmall) then
!
! integrals ... inlined ... for speed ...
!
           ACTOT = ZERO
           y = ONE/SPSIG
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             do is = 2, msc
               ACTOT = ACTOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
           end do
           !tmp = ONE/SPSIG
           !ACTOT       = DINTSPEC_Y(IP,ACLOC,tmp)
           ETOT_SPSIG = ZERO
           y = SIGPOW(:,1) 
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             do is = 2, msc
               ETOT_SPSIG = ETOT_SPSIG + ONEHALF*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
           end do
           !tmp = SIGPOW(:,1)
           !ETOT_SPSIG  = DINTSPEC_Y(IP,ACLOC,tmp)
           ETOT_WK = ZERO
           y = WK(:,IP) 
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             do is = 2, msc
               ETOT_WK = ETOT_WK + ONEHALF*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
           end do
           !tmp = WK(:,IP)
           !ETOT_WK     = DINTSPEC_Y(IP,ACLOC,tmp)
           ETOT_ISQ_WK = ZERO 
           y = ONE/SQRT(WK(:,IP))
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             do is = 2, msc
               ETOT_ISQ_WK = ETOT_ISQ_WK + ONEHALF*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
           end do
           !tmp = ONE/SQRT(WK(:,IP))
           !ETOT_ISQ_WK = DINTSPEC_Y(IP,ACLOC,tmp)
           ETOT_SQ_WK = ZERO
           y = SQRT(WK(:,IP)) 
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             do is = 2, msc
               ETOT_SQ_WK = ETOT_SQ_WK + ONEHALF*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
           end do
           !tmp = SQRT(WK(:,IP))
           !ETOT_SQ_WK  = DINTSPEC_Y(IP,ACLOC,tmp) 
!
! tail factors ...
!
           DS          = SPSIG(MSC) - SPSIG(MSC-1)

           ATAIL       = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,1) * DDIR * DS
           ETAIL       = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,2) * DDIR * DS
           ESIGTAIL    = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,3) * DDIR * DS
!
! tail factors ... borowed from SWAN
!
           ACTOT       = ACTOT        + PTAIL(5)  * ATAIL 
           ETOT        = ETOT         + PTAIL(6)  * ETAIL 
           ETOT_SPSIG  = ETOT_SPSIG   + PTAIL(7)  * ETAIL 
           ETOT_ISQ_WK = ETOT_ISQ_WK  + PTAIL(5)  * ETAIL / (SQRT(WK(MSC,IP)))
           ETOT_SQ_WK  = ETOT_SQ_WK   + PTAIL(5)  * ETAIL * (SQRT(WK(MSC,IP)))
           ETOT_WK     = ETOT_WK      + PTAIL(8)  * ETAIL * WK(MSC,IP)
!
! integral parameters ...
!
           KME01       =  ETOT_WK / ETOT
           KMWAM       = (ETOT/ETOT_ISQ_WK)**2
           KMWAM2      = (ETOT_SQ_WK/ETOT)**2

         else
!
! no or too less energy ...
!
           KME01       = 10.0_rkind
           KMWAM       = 10.0_rkind
           KMWAM2      = 10.0_rkind

         end if 
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
