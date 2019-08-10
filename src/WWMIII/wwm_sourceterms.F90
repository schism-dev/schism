#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
!2do add mean quantities for 
      SUBROUTINE SOURCETERMS (IP, ACLOC, IMATRA, IMATDA, LRECALC, ISELECT, CALLFROM)
      
      USE DATAPOOL, ONLY : MSC, MDC, MNP, WK, LINID, THR, UFRIC
      USE DATAPOOL, ONLY : CD, TAUTOT, TAUWX, TAUWY, DEP, AC1, AC2
      USE DATAPOOL, ONLY : PI2, CG, G9, ZERO, ALPHA_CH, QBLOCAL
      USE DATAPOOL, ONLY : USTDIR, Z0, SMALL, VERYSMALL, MSC_HF
      USE DATAPOOL, ONLY : DDIR, SPSIG, SPDIR, FRINTF, LMAXETOT, LADVTEST
      USE DATAPOOL, ONLY : USTDIR, Z0, SMALL, VERYSMALL, MSC_HF, DDIR
      USE DATAPOOL, ONLY : SPSIG, SPDIR, FRINTF, ONE, RHOA, RHOAW, TAUHF
      USE DATAPOOL, ONLY : TAUW, FR, MESNL, MESIN, MESDS, MESBF, MESBR
      USE DATAPOOL, ONLY : MESTR, ISHALLOW, DS_INCR, IOBP, IOBPD, MEVEG
      USE DATAPOOL, ONLY : LNANINFCHK, DBG, IFRIC, RTIME, DISSIPATION
      USE DATAPOOL, ONLY : AIRMOMENTUM, ONEHALF, NSPEC, RKIND, ICOMP
#ifdef WWM_MPI
      USE DATAPOOL, ONLY : myrank
#endif
         USE SdsBabanin
#ifdef SNL4_TSA
         USE W3SNLXMD
#endif

#ifdef ST41
         USE W3SRC4MD_OLD
#elif ST42
         USE W3SRC4MD
#endif
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP, ISELECT

         REAL(rkind), INTENT(OUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind), INTENT(IN)  :: ACLOC(MSC,MDC)

         LOGICAL, INTENT(IN) :: LRECALC
         CHARACTER(LEN=*), INTENT(IN) :: CALLFROM

         INTEGER        :: IS, ID, IS0, IK, ITH, IDISP, JU, NZZ
         REAL(rkind)    :: WIND10, WINDTH
         REAL(rkind)    :: FPM
         REAL(rkind)    :: SME01, SME10, KME01, KMWAM, KMWAM2
         REAL(rkind)    :: SME01WS, SME10WS
         REAL(rkind)    :: HS, ETOT, FPMH,FPM4 
         REAL(rkind)    :: LPOW, MPOW, a1, a2, ETOTWS
         REAL(rkind)    :: XRR, XPP, XFLT, XREL, FACP, XFILT
         REAL(rkind)    :: TEMP2(MSC), SSBRL(MSC,MDC)
         REAL(rkind)    :: AWW3(NSPEC), IMATRAWW3(MSC,MDC), IMATDAWW3(MSC,MDC)
         REAL(rkind)    :: EDENS(MSC), KDS(MSC), ABAB(MSC)
         REAL(rkind)    :: WHITECAP(1:4),AKMEAN,XKMEAN,F1MEAN,TMPAC(MDC,MSC),TEMP(MSC), FCONST(MSC)
         REAL(rkind)    :: ACLOC1(MSC,MDC)


         REAL(rkind)    :: SSNL3(MSC,MDC), SSNL4(MSC,MDC), SSINL(MSC,MDC), SSDS(MSC,MDC), DSSNL4(MSC,MDC)
         REAL(rkind)    :: SSVEG(MSC,MDC),DSSVEG(MSC,MDC)
         REAL(rkind)    :: SSBF(MSC,MDC), SSBR(MSC,MDC), SSINE(MSC,MDC), DSSNL3(MSC,MDC), DSSBR(MSC,MDC)
         REAL(rkind)    :: TMP_IN(MSC), TMP_DS(MSC), WN2(MSC*MDC),SPRDD(MDC),AK2VGM1,AKM1

         REAL(rkind)    :: EMEAN, FMEAN, FMEAN1, WNMEAN, AMAX, AS
         REAL(rkind)    :: FMEANWS, TAUWAX, TAUWAY, XJ, FLLOWEST, GADIAG
         REAL(rkind)    :: IMATDA1D(NSPEC), IMATRA1D(NSPEC), SUMACLOC, IMATRAT(MSC,MDC), BRLAMBDA(NSPEC)
         REAL(rkind)    :: IMATRA_WAM(MDC,MSC), IMATDA_WAM(MDC,MSC), TAILFACTOR, FLOGSPRDM1, SNL3(MSC,MDC), DSNL3(MSC,MDC)

#ifdef SNL4_TSA
         REAL           :: IMATRA_TSA(MDC,MSC), IMATDA_TSA(MDC,MSC), TMPAC_TSA(MDC,MSC), CG_TSA(MSC), WK_TSA(MSC), DEP_TSA
#endif

         REAL           :: XNL(MSC,MDC), DDIAG(MSC,MDC), ACLOC_WRT(MSC,MDC), DEP_WRT, SPSIG_WRT(MSC), SPDIR_WRT(MDC)
         INTEGER        :: IERR_WRT

#ifdef TIMINGS 
         REAL(rkind)        :: T1, T2
         REAL(rkind), SAVE  :: TIME1, TIME2, TIME3, TIME4, TIME5, TIME6, TIME7, TIME8, TIME9
#endif

         LOGICAL        :: LWINDSEA(MSC,MDC)
         REAL(rkind)    :: XLCKS(MDC,MSC)

#ifdef TIMINGS
         call WAV_MY_WTIME(TIME1)
#endif 

!         IF (LMAXETOT .AND. .NOT. LADVTEST .AND. ISHALLOW(IP) .EQ. 1 .AND. .NOT. LRECALC) THEN
!           CALL BREAK_LIMIT(IP,ACLOC,SSBRL) ! Miche to reduce stiffness of source terms ...
!         END IF

         IDISP = 999

         IF (.NOT. LRECALC) THEN
           ACLOC1 = AC1(:,:,IP)
           CALL MEAN_WAVE_PARAMETER(IP,ACLOC1,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2) ! 1st guess ... 
         END IF

         SUMACLOC    = SUM(ACLOC)
!
         WIND10      = ZERO
         BRLAMBDA    = ZERO
         SSINE       = zero
         SSINL       = zero
         SSNL3       = zero
         SSNL4       = zero
         DSSNL4      = zero
         SSBR        = zero
         SSBF        = zero
         SSVEG       = zero
         DSSVEG      = zero
         IMATRA_WAM  = zero
         IMATDA_WAM  = zero
         TMPAC       = zero
         IMATRA      = zero 
         IMATDA      = zero 
         IMATRAWW3   = zero 
         IMATDAWW3   = zero 
         QBLOCAL(IP) = zero 

#ifdef ST_DEF
         IF (MESDS == 1 .OR. MESIN .EQ. 1) THEN
           DO IS = 1, MSC
             DO ID = 1, MDC
               AWW3(ID + (IS-1) * MDC) = ACLOC(IS,ID) * CG(IS,IP)
             END DO
           END DO
         END IF
         XPP     = 0.15
         XRR     = 0.10
         XFILT  = 0.05
         XPP     = MAX ( 1.E-6_rkind , XPP )
         XRR     = MAX ( 1.E-6_rkind , XRR )
         XREL   = XRR
         XFILT  = MAX ( ZERO , XFILT )
         XFLT   = XFILT
         FACP   = 2*XPP / PI2 * 0.62E-3 * PI2**4 / G9**2
         DO IK=1, NK
           WN2(1+(IK-1)*NTH) = WK(IK,IP)
         END DO
         DO IK=1, NK
           IS0    = (IK-1)*NTH
           DO ITH=2, NTH
             WN2(ITH+IS0) = WN2(1+IS0)
           END DO
         END DO
#endif

#ifdef DEBUG
         WRITE(12,*) '0', SUM(IMATRA), SUM(IMATDA), ISELECT, MESIN, CALLFROM
#endif

#ifdef TIMINGS
         call WAV_MY_WTIME(TIME2)
#endif 
         IF ((ISELECT .EQ. 1 .OR. ISELECT .EQ. 10 .OR. ISELECT .EQ. 20) .AND. .NOT. LRECALC) THEN
           IF (IOBP(IP) .EQ. 0) THEN
             IF (MESIN == 1) THEN ! Ardhuin et al. 2010
               CALL SET_WIND( IP, WIND10, WINDTH )
               TAUWX(IP) = ZERO
               TAUWY(IP) = ZERO               
               IF (SUMACLOC .LT. THR .AND. WIND10 .GT. THR) THEN
                 CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
                 CALL SIN_LIN_CAV(IP,WIND10,WINDTH,FPM,IMATRA,SSINL)
                 !WRITE(*,'(I10,4F20.10)') IP, WIND10, FPM, SUM(IMATRA)
               ELSE
                 MSC_HF(IP) = MSC
                 AS      = 0. 
#ifdef ST41
                 CALL W3SPR4_OLD ( AWW3, CG(:,IP), WK(:,IP), EMEAN, FMEAN, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
                 CALL W3SIN4_OLD ( AWW3, CG(:,IP), WK(:,IP), WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, IMATRA1D, IMATDA1D, LLWS)
                 CALL W3SPR4_OLD ( AWW3, CG(:,IP), WK(:,IP), EMEAN, FMEAN, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
                 CALL W3SIN4_OLD ( AWW3, CG(:,IP), WK(:,IP), WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, IMATRA1D, IMATDA1D, LLWS)
#elif ST42
                 CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), EMEAN, FMEAN, FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
                 CALL W3SIN4 ( IP, AWW3, CG(:,IP), WN2,  WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, IMATRA1D, IMATDA1D, LLWS, BRLAMBDA)
                 CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), EMEAN, FMEAN, FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)  
                 CALL W3SIN4 ( IP, AWW3, CG(:,IP), WN2, WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, IMATRA1D, IMATDA1D, LLWS, BRLAMBDA) 
#else
                 WRITE(DBG%FHNDL,*) 'NO ST42 or ST41 chosen but MESIN == 1'
                 CALL WWM_ABORT('stop wwm_sourceterms l.176')
#endif
                 CALL ONED2TWOD(IMATRA1D,IMATRAWW3)
                 SSINE = IMATRAWW3
                 CALL ONED2TWOD(IMATDA1D,IMATDAWW3)
                 DO ID = 1, MDC
                   IMATRA(:,ID) = IMATRAWW3(:,ID) / CG(:,IP) + IMATRA(:,ID)
                   IMATDA(:,ID) = ZERO!IMATDAWW3(:,ID) !/ CG(:,IP) 
                 END DO
               END IF
               IF (LNANINFCHK) THEN
                 IF (SUM(IMATRA) .NE. SUM(IMATRA) .OR. SUM(IMATDA) .NE. SUM(IMATDA)) THEN
                   WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT NORMAL', IP, 'DUE TO SIN', SUM(IMATRA), SUM(IMATDA)
                   CALL WWM_ABORT('NAN AT wwm_sourceterms.F90 l.189')
                 END IF
               ENDIF
             ELSE IF (MESIN == 2) THEN ! Cycle 4, Bidlot et al. ...
               CALL WWM_ABORT('PLEASE USE "LSOURCEWAM = T" FOR ECWAM SOURCE TERM FORMULATION') 
             ELSE IF (MESIN == 3) THEN ! Makin & Stam
               CALL SET_WIND( IP, WIND10, WINDTH )
               IFRIC = 4
               CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
               CALL SIN_LIN_CAV( IP, WIND10, WINDTH, FPM, IMATRA, SSINL)
               CALL SIN_MAKIN( IP, WIND10, WINDTH, KME01,ETOT,ACLOC,IMATRA,IMATDA,SSINE)
             ELSE IF (MESIN == 4) THEN ! Donealan et al.
               CALL SET_WIND( IP, WIND10, WINDTH )
               IFRIC = 1
               CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
               IF (.NOT. LINID) CALL SIN_LIN_CAV( IP, WIND10, WINDTH,FPM,IMATRA,SSINL)
               CALL SWIND_DBYB (IP,WIND10,WINDTH,IMATRA,SSINE)
             ELSE IF (MESIN == 5) THEN ! Cycle 3
               CALL SET_WIND( IP, WIND10, WINDTH )
               IFRIC = 1
               CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
               IF (.NOT. LINID) CALL SIN_LIN_CAV( IP, WIND10, WINDTH,FPM,IMATRA,SSINL)
               CALL SIN_EXP_KOMEN( IP, WINDTH, ACLOC, IMATRA, IMATDA, SSINE )
             END IF ! MESIN
           END IF ! IOBP
           !IF (IDISP == IP) THEN
             !CALL PLOT_SHADED_CONTOUR_POLAR(SPSIG/PI2,SPDIR*RADDEG,MSC,MDC,IMATRA/MAXVAL(IMATRA),50,MSC,MDC,'BEFORE ANY CALL')
             !PAUSE
           !END IF
         END IF ! ISELECT 

#ifdef DEBUG
         WRITE(12,*) '1', SUM(IMATRA), SUM(IMATDA)
#endif

#ifdef TIMINGS
         call WAV_MY_WTIME(TIME3)
#endif 
         IF (LNANINFCHK) THEN
           IF (SUM(IMATRA) .NE. SUM(IMATRA) .OR. SUM(IMATDA) .NE. SUM(IMATDA)) THEN
             WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, '   DUE TO SIN', SUM(IMATRA), SUM(IMATDA)
             CALL WWM_ABORT('NAN in wwm_sourceterm.F90 l.311')
           END IF
         ENDIF

         IF ((ISELECT.EQ.2 .OR. ISELECT.EQ.10 .OR. ISELECT.EQ.20) .AND. .NOT. LRECALC) THEN
           IF (IOBP(IP) .EQ. 0) THEN
             IF (MESNL .EQ. 1) THEN
               CALL SNL41 (IP, KMWAM, ACLOC, IMATRA, IMATDA, SSNL4, DSSNL4)
             ELSE IF (MESNL .EQ. 2) THEN
               CALL SNL4(IP, KMWAM, ACLOC, IMATRA, IMATDA)
             ELSE IF (MESNL .EQ. 3) THEN
               CALL SNL42(IP, KMWAM, ACLOC, IMATRA, IMATDA)
             ELSE IF (MESNL .EQ. 4) THEN
               CALL SNL43(IP, KMWAM, ACLOC, IMATRA, IMATDA)
             ELSE IF (MESNL .EQ. 5) THEN
               ACLOC_WRT = REAL(ACLOC)
               SPSIG_WRT = REAL(SPSIG)
               SPDIR_WRT = REAL(SPDIR)
               DEP_WRT   = DEP(IP)
               CALL WWMQUAD_WRT (ACLOC_WRT,SPSIG_WRT,SPDIR_WRT,MDC,MSC,DEP_WRT,1,XNL,DDIAG,IERR_WRT)
               IF (IERR_WRT .GT. 0) THEN
                 WRITE (DBG%FHNDL,*) 'XNL_WRT ERROR', IERR_WRT
                 CALL WWM_ABORT('XNL_WRT ERROR')
               ELSE
                 IMATRA(:,:) = IMATRA(:,:) + XNL (:,:)
                 IMATDA(:,:) = IMATDA(:,:) + DDIAG(:,:)
               END IF
             ELSE IF (MESNL .EQ. 6) THEN
#ifdef SNL4_TSA
               DO IS = 1, MSC
                 DO ID = 1, MDC
                   TMPAC_TSA(ID,IS) = ACLOC(IS,ID) * CG(IS,IP)
                 END DO
               END DO
               CG_TSA = CG(:,IP)
               WK_TSA = WK(:,IP)
               DEP_TSA = DEP(IP)
               NZZ = (MSC*(MSC+1))/2
               CALL W3SNLX ( TMPAC_TSA, CG_TSA, WK_TSA, DEP_TSA, NZZ, IMATRA_TSA, IMATDA_TSA)
               DO IS = 1, MSC
                 DO ID = 1, MDC
                   IMATRA(IS,ID) = IMATRA(IS,ID) + IMATRA_TSA(ID,IS) / CG(IS,IP)
                   IMATDA(IS,ID) = IMATDA(IS,ID) + IMATDA_TSA(ID,IS)
                 END DO
               END DO
#endif
             END IF
           END IF
           IF (IDISP == IP) THEN
             !CALL PLOT_SHADED_CONTOUR_POLAR(SPSIG/PI2,SPDIR*RADDEG,MSC,MDC,IMATRA/MAXVAL(IMATRA),50,MSC,MDC,'BEFORE ANY CALL')
             !PAUSE
           END IF
         END IF ! ISELECT

         IF (LNANINFCHK) THEN
           IF (SUM(IMATRA) .NE. SUM(IMATRA) .OR. SUM(IMATDA) .NE. SUM(IMATDA)) THEN
             WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, 'DUE TO SNL4'
             CALL WWM_ABORT('NAN at wwm_sourceterms.F90 l.368')
           END IF
         ENDIF

#ifdef DEBUG
         WRITE(12,*) '2', SUM(IMATRA), SUM(IMATDA)
#endif


#ifdef TIMINGS
         call WAV_MY_WTIME(TIME4)
#endif 
         IF ((ISELECT.EQ.3 .OR. ISELECT.EQ.10 .OR. ISELECT.EQ.20) .AND. .NOT. LRECALC) THEN

           IMATRAT = IMATRA

           IF (IOBP(IP) .EQ. 0 .OR. IOBP(IP) .EQ. 4) THEN

             IF (MESDS == 1) THEN
#ifdef ST41
               CALL W3SDS4_OLD(AWW3,WK(:,IP),CG(:,IP),UFRIC(IP),USTDIR(IP),DEP(IP),IMATRA1D, IMATDA1D) 
#elif ST42
               CALL W3SDS4(AWW3,WK(:,IP),CG(:,IP),UFRIC(IP),USTDIR(IP),DEP(IP),IMATRA1D,IMATDA1D,BRLAMBDA,WHITECAP)
#endif
               CALL ONED2TWOD(IMATRA1D,IMATRAWW3)
               CALL ONED2TWOD(IMATDA1D,IMATDAWW3)
               DO ID = 1, MDC
                 SSDS(:,ID)   = IMATRAWW3(:,ID) / CG(:,IP)
                 IMATRA(:,ID) = IMATRA(:,ID)+IMATRAWW3(:,ID) / CG(:,IP)
                 !IMATDA(:,ID) = IMATDA(:,ID)+IMATDAWW3(:,ID) !/ CG(:,IP)
               END DO
             ELSE IF (MESDS == 2) THEN
               CALL WWM_ABORT('PLEASE USE LSOURCEWAM FOR ECWAM SOURCE TERM FORMULATION')
             ELSE IF (MESDS == 3) THEN
               CALL SDS_NEDWAM_CYCLE4( IP, KMWAM, SME01, ETOT, ACLOC, IMATRA, IMATDA, SSDS  ) ! NEDWAM 
             ELSE IF (MESDS == 4) THEN
               ABAB = 1.
               LPOW = 4.
               MPOW = 4.
               a1  = 0.00000045
               a2  = 0.0000095
!              LPOW = 2.
!              MPOW = 2.
!              a1  = 0.0002
!              a2  = 0.003
!  0.0002 0.003 2.0 2.0 KOM
!  0.00000045 0.0000095 4.0 4.0 BD
               DO IS = 1, MSC
                 EDENS(IS) = 0.
                 DO ID = 1, MDC
                   EDENS(IS) = EDENS(IS) + ACLOC(IS,ID) *  SPSIG(IS) * PI2 * DDIR
                 END DO
               END DO
               CALL CALC_SDS(IP,MSC,EDENS,FR,Kds,ABAB,LPOW,MPOW,a1,a2,ACLOC,IMATRA,IMATDA,SSDS)
!              CALL SSWELL(IP,ETOT,ACLOC,IMATRA,IMATDA,URMSTOP,CG0)
             ELSE IF (MESDS == 5) THEN
               CALL SDS_CYCLE3 ( IP, KMWAM, SME10, ETOT, ACLOC, IMATRA, IMATDA, SSDS ) 
             END IF
           END IF
#ifdef VDISLIN
           IF (IDISP == IP) THEN
             CALL PLOT_SHADED_CONTOUR_POLAR(SPSIG/PI2,SPDIR*RADDEG,MSC,MDC,IMATRAT-IMATRA,50,MSC,MDC,'BEFORE ANY CALL')
           END IF
#endif
         END IF

#ifdef DEBUG
         WRITE(12,*) '3', SUM(IMATRA), SUM(IMATDA)
#endif

#ifdef TIMINGS
         call WAV_MY_WTIME(TIME5)
#endif 
         IF (((ISELECT.EQ.4 .OR. ISELECT.EQ.10) .AND.ISHALLOW(IP).EQ.1) .AND. .NOT. LRECALC) THEN
           IF (SUMACLOC .GT. VERYSMALL) THEN
             IF ((MESTR .EQ. 1).or.(MESTR .eq. 5) ) THEN
               CALL TRIAD_ELDEBERKY(IP,HS,SME01,ACLOC,IMATRA,IMATDA,SSNL3,DSSNL3)
             ELSE IF (MESTR .EQ. 2) THEN
               CALL SNL31 (IP,HS,SME01,ACLOC,IMATRA,IMATDA,SSNL3)
             ELSE IF (MESTR .EQ. 3) THEN
               CALL SNL32 (IP,HS,SME01,ACLOC,IMATRA,IMATDA,SSNL3)
             ELSE IF (MESTR .EQ. 4) THEN
               CALL SNL33 (IP,HS,SME01,ACLOC,IMATRA,IMATDA,SSNL3)
             ELSE IF (MESTR .EQ. 6) THEN
               CALL WWM_ABORT('NOT READ YET')
               CALL TRIAD_DINGEMANS (IP,ACLOC,IMATRA,IMATDA,SSNL3)
             ELSE IF (MESTR .EQ. 7) THEN
               CALL WWM_ABORT('NOT READ YET')
               CALL snl3ta(ip,snl3,dsnl3)
               SSNL3 = SNL3
               IMATRA = IMATRA + SNL3
               IMATDA = IMATDA + DSNL3
             END IF
           END IF
         END IF

#ifdef DEBUG
        WRITE(12,*) '4', SUM(IMATRA), SUM(IMATDA)
#endif


#ifdef TIMINGS
         call WAV_MY_WTIME(TIME6)
#endif 
         IF (MESBR .EQ. 1) THEN
           IF (((ISELECT.EQ.5 .OR. ISELECT.EQ.10 .OR. ISELECT.EQ.30) .AND. ISHALLOW(IP) .EQ. 1) .AND. .NOT. LRECALC) THEN
             IF (IOBP(IP) == 0 .OR. IOBP(IP) == 4 .OR. IOBP(IP) == 3) THEN
               CALL SDS_SWB(IP,SME01,KMWAM,ETOT,HS,ACLOC,IMATRA,IMATDA,SSBR,DSSBR)
             END IF
           ENDIF
         END IF
#ifdef TIMINGS
         call WAV_MY_WTIME(TIME7)
#endif 

#ifdef DEBUG
         WRITE(12,*) '5', SUM(IMATRA), SUM(IMATDA)
#endif

         IF (MESBF .GE. 1) THEN
           IF (((ISELECT.EQ.6 .OR. ISELECT.EQ.10 .OR. ISELECT.EQ.30) .AND. ISHALLOW(IP) .EQ. 1) .AND. .NOT. LRECALC) THEN
             IF (IOBP(IP) == 0 .OR. IOBP(IP) == 4 .OR. IOBP(IP) == 3) THEN
              CALL SDS_BOTF(IP,ACLOC,IMATRA,IMATDA,SSBR,DSSBR)
             END IF
           END IF
         ENDIF

         IF (MEVEG .GE. 1) THEN
           IF (((ISELECT.EQ.6 .OR. ISELECT.EQ.10 .OR. ISELECT.EQ.30) .AND. ISHALLOW(IP) .EQ. 1) .AND. .NOT. LRECALC) THEN
             IF (IOBP(IP) == 0 .OR. IOBP(IP) == 4 .OR. IOBP(IP) == 3) THEN
              CALL VEGDISSIP (IP,IMATRA,IMATDA,SSVEG,DSSVEG,ACLOC,DEP(IP),ETOT,SME01,KME01)
             END IF
           END IF
         ENDIF

         IF (LNANINFCHK) THEN
           IF (SUM(IMATRA) .NE. SUM(IMATRA) .OR. SUM(IMATDA) .NE. SUM(IMATDA)) THEN
             WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, '   DUE TO SBF' 
             CALL WWM_ABORT('NAN in wwm_sourceterms.F90 at l.419')
           END IF
         ENDIF

#ifdef DEBUG
         WRITE(12,*) '6', SUM(IMATRA), SUM(IMATDA)
#endif

#ifdef TIMINGS
        call WAV_MY_WTIME(TIME8)
#endif 
!------------------------------------------------------------------------------------------------------------------------!
!-------------------------------- RECALCULATE ALL SOURCE TERMS BASED ON THE NEW SPECTRA ---------------------------------! 
!------------------------------------------------------------------------------------------------------------------------!
         IF (LRECALC .and. IOBP(IP) .EQ. 0) THEN

           IF (MESIN == 1) THEN
             AS      = 0.
#ifdef ST41
             CALL W3SPR4_OLD ( AWW3, CG(:,IP), WK(:,IP), EMEAN, FMEAN, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
             CALL W3SIN4_OLD ( AWW3, CG(:,IP), WK(:,IP), WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, IMATRA1D, IMATDA1D, LLWS)
             CALL W3SPR4_OLD ( AWW3, CG(:,IP), WK(:,IP), EMEAN, FMEAN, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
#elif ST42
             CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), EMEAN, FMEAN, FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
             CALL W3SIN4 ( IP, AWW3, CG(:,IP), WN2,  WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, IMATRA1D, IMATDA1D, LLWS, BRLAMBDA)
             CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), EMEAN, FMEAN, FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
#else
             WRITE(DBG%FHNDL,*) 'NO ST42 or ST41 chosen but MESIN == 1'
             CALL WWM_ABORT('stop wwm_sourceterms l.186')
#endif
           ELSEIF (MESIN == 2) THEN
           ELSEIF (MESIN == 3) THEN
           ELSEIF (MESIN == 4) THEN
           ELSEIF (MESIN == 5) THEN
           ENDIF

           DISSIPATION(IP) = 0.
           AIRMOMENTUM(IP) = 0.
           DO ID = 1, MDC
             TMP_DS = ( SSBR(:,ID) + SSBF(:,ID) + SSDS(:,ID) ) * SPSIG * DDIR
             TMP_IN = ( SSINE(:,ID) + SSINL(:,ID) ) * SPSIG * DDIR
             DO IS = 2, MSC
               DISSIPATION(IP) = DISSIPATION(IP) + ONEHALF * ( TMP_DS(IS) + TMP_DS(IS-1) ) * DS_INCR(IS)
               AIRMOMENTUM(IP) = AIRMOMENTUM(IP) + ONEHALF * ( TMP_IN(IS) + TMP_IN(IS-1) ) * DS_INCR(IS)
             END DO
           END DO

         ENDIF

#ifdef TIMINGS 
         call WAV_MY_WTIME(TIME9)
#ifdef WWM_MPI 
         IF (IP == MNP .AND. myrank == 0 ) THEN
#else 
         IF (IP == MNP) THEN
#endif
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----SOURCE TIMINGS-----'
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'PREPARATION        ', TIME2-TIME1
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SIN                ', TIME3-TIME2
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SDS                ', TIME4-TIME3
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SNL4               ', TIME5-TIME4
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SNL3               ', TIME6-TIME5
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SBR                ', TIME7-TIME6
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SBF                ', TIME8-TIME7
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'RECALC             ', TIME9-TIME8
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'TOTAL              ', TIME9-TIME1
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '------END-TIMINGS-  ---'
         ENDIF
#endif 

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************

