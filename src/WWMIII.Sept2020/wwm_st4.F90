#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ST4_PRE (IP, WALOC, PHI, DPHIDN, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
      USE DATAPOOL
      USE W3SRC4MD
      IMPLICIT NONE

      INTEGER, INTENT(IN)        :: IP
      REAL(rkind), INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)

      REAL(rkind), INTENT(OUT)   :: PHI(NUMSIG,NUMDIR), DPHIDN(NUMSIG,NUMDIR)
      REAL(rkind), INTENT(OUT)   :: SSINE(NUMSIG,NUMDIR), DSSINE(NUMSIG,NUMDIR) 
      REAL(rkind), INTENT(OUT)   :: SSDS(NUMSIG,NUMDIR),DSSDS(NUMSIG,NUMDIR)
      REAL(rkind), INTENT(OUT)   :: SSNL4(NUMSIG,NUMDIR),DSSNL4(NUMSIG,NUMDIR)
      REAL(rkind), INTENT(OUT)   :: SSINL(NUMSIG,NUMDIR)

      INTEGER      :: IS, ID, ITH, IK, IS0

      REAL(rkind)  :: AWW3(NSPEC)
      REAL(rkind)  :: VDDS(NSPEC), VSDS(NSPEC), BRLAMBDA(NSPEC)
      REAL(rkind)  :: WN2(NUMSIG*NUMDIR), WHITECAP(1:4)
      REAL(rkind)  :: VSIN(NSPEC), VDIN(NSPEC)

      REAL(rkind)  :: ETOT, FAVG, FMEAN1, WNMEAN, AS, SUMWALOC, FAVGWS
      REAL(rkind)  :: TAUWAX, TAUWAY, AMAX, FPM, WIND10, WINDTH
      REAL(rkind)  :: HS,SME01,SME10,KME01,KMWAM,KMWAM2

      DO IS = 1, NUMSIG
        DO ID = 1, NUMDIR
          AWW3(ID + (IS-1) * NUMDIR) = WALOC(IS,ID) * CG(IS,IP)
        END DO
      END DO

      DO IK=1, NK
        WN2(1+(IK-1)*NTH) = WK(IK,IP)
      END DO

      DO IK=1, NK
        IS0    = (IK-1)*NTH
        DO ITH=2, NTH
          WN2(ITH+IS0) = WN2(1+IS0)
        END DO
      END DO
!
! wind input
!
      TAUWX(IP)  = ZERO
      TAUWY(IP)  = ZERO               
      SSINL      = ZERO
      NUMSIG_HF(IP) = NUMSIG
      AS         = 0.
      BRLAMBDA   = ZERO

      IF (MESIN .GT. 0) THEN

        CALL SET_WIND( IP, WIND10, WINDTH )
        CALL SET_FRICTION( IP, WALOC, WIND10, WINDTH, FPM )
        LLWS=.TRUE.
#ifdef DEBUGSRC
        WRITE(740+myrank,*) '1: input value USTAR=', UFRIC(IP), ' USTDIR=', USTDIR(IP)
#endif
        CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), EMEAN(IP), FMEAN(IP), FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS(IP))
#ifdef DEBUGSRC
        WRITE(740+myrank,*) '1: out value USTAR=', UFRIC(IP), ' USTDIR=', USTDIR(IP)
        WRITE(740+myrank,*) '1: out value EMEAN=', EMEAN(IP), ' FMEAN=', FMEAN(IP)
        WRITE(740+myrank,*) '1: out value FMEAN1=', FMEAN1, ' WNMEAN=', WNMEAN
        WRITE(740+myrank,*) '1: out value CD=', CD(IP), ' Z0=', Z0(IP)
        WRITE(740+myrank,*) '1: out value ALPHA=', ALPHA_CH(IP), ' FMEANWS=', FMEANWS(IP)
#endif
        IF (EMEAN(IP) .LT. THR .AND. WIND10 .GT. THR) CALL SIN_LIN_CAV(IP,WIND10,WINDTH,FPM,SSINL)
        CALL W3SIN4 ( IP, AWW3, CG(:,IP), WN2, WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, VSIN, VDIN, LLWS, BRLAMBDA)
#ifdef DEBUGSRC
        WRITE(740+myrank,*) '1: WINDTH=', WINDTH, ' Z0=', Z0(IP), ' CD=', CD(IP)
        WRITE(740+myrank,*) '1: UFRIC=', UFRIC(IP), 'WIND10=', WIND10, ' RHOAW=', RHOAW
        WRITE(740+myrank,*) '1: TAUWX=', TAUWX(IP), ' TAUWY=', TAUWY(IP)
        WRITE(740+myrank,*) '1: TAUWAX=', TAUWAX, ' TAUWAY=', TAUWAY
        WRITE(740+myrank,*) '1: W3SIN4min/max/sum(VSIN)=', minval(VSIN), maxval(VSIN), sum(VSIN)
        WRITE(740+myrank,*) '1: W3SIN4min/max/sum(VDIN)=', minval(VDIN), maxval(VDIN), sum(VDIN)
#endif
        CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), EMEAN(IP), FMEAN(IP), FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS(IP))  
        CALL W3SIN4 ( IP, AWW3, CG(:,IP), WN2, WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, VSIN, VDIN, LLWS, BRLAMBDA)
#ifdef DEBUGSRC
        WRITE(740+myrank,*) '2: W3SIN4min/max/sum(VSIN)=', minval(VSIN), maxval(VSIN), sum(VSIN)
        WRITE(740+myrank,*) '2: W3SIN4min/max/sum(VDIN)=', minval(VDIN), maxval(VDIN), sum(VDIN)
#endif
        CALL CONVERT_VS_VD_WWM(IP, VSIN, VDIN, SSINE, DSSINE)
      ENDIF

      IF (MESNL .GT. 0) THEN
         CALL MEAN_WAVE_PARAMETER(IP,WALOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)
         CALL DIASNL4WW3(IP, KMWAM, WALOC, SSNL4, DSSNL4)
      END IF

      IF (MESDS .GT. 0) THEN
        CALL W3SDS4(AWW3,WK(:,IP),CG(:,IP),UFRIC(IP),USTDIR(IP),DEP(IP),VSDS,VDDS,BRLAMBDA,WHITECAP)

#ifdef DEBUGSRC
        WRITE(740+myrank,*) '2: W3SDS4min/max/sum(VSDS)=', minval(VSDS), maxval(VSDS), sum(VSDS)
        WRITE(740+myrank,*) '2: W3SDS4min/max/sum(VDDS)=', minval(VDDS), maxval(VDDS), sum(VDDS)
#endif
        CALL CONVERT_VS_VD_WWM(IP, VSDS, VDDS, SSDS, DSSDS)
      ENDIF
!
      PHI    = SSINL + SSINE  + SSNL4  + SSDS
      DPHIDN =         DSSINE + DSSNL4 + DSSDS
!
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ST4_POST (IP, WALOCOLD, WALOC)

        USE DATAPOOL
        USE W3SRC4MD
        IMPLICIT NONE
       
        INTEGER, INTENT(IN)        :: IP
        REAL(rkind), INTENT(IN)    :: WALOCOLD(NUMSIG,NUMDIR)
        REAL(rkind), INTENT(OUT)   :: WALOC(NUMSIG,NUMDIR)
       
        REAL(rkind)                :: SSINE(NUMSIG,NUMDIR),DSSINE(NUMSIG,NUMDIR)
        REAL(rkind)                :: SSDS(NUMSIG,NUMDIR),DSSDS(NUMSIG,NUMDIR)
        REAL(rkind)                :: SSINL(NUMSIG,NUMDIR)

        INTEGER                    :: IS, ID, IK, ITH, ITH2, IS0, NKH, NKH1

        REAL(rkind)                :: PHI(NUMSIG,NUMDIR), DPHIDN(NUMSIG,NUMDIR)
        REAL(rkind)                :: AWW3(NSPEC), AWW3OLD(NSPEC), WN2(NUMSIG*NUMDIR), BRLAMBDA(NSPEC)
        REAL(rkind)                :: SSINE_WW3(NSPEC), DSSINE_WW3(NSPEC), TMP_DS(NUMSIG)
        REAL(rkind)                :: SSDS_WW3(NSPEC), DSSDS_WW3(NSPEC), SSNL4_WW3(NSPEC), DSSNL4_WW3(NSPEC)
        REAL(rkind)                :: SSBF_WW3(NSPEC), DSSBF_WW3(NSPEC), SSNL4(NUMSIG,NUMDIR), DSSNL4(NUMSIG,NUMDIR)
        REAL(rkind)                :: SSBF(NUMSIG,NUMDIR), DSSBF(NUMSIG,NUMDIR)

        REAL(rkind)                :: ETOT, FAVG, FMEAN1, WNMEAN, AS, FAVGWS
        REAL(rkind)                :: TAUWAX, TAUWAY, AMAX, WIND10, WINDTH, EBAND, B1BAND
        REAL(rkind)                :: WHITECAP(1:4), SUMWALOC, FPM, FH1, FH2, TAUBBL(2)
        REAL(rkind)                :: PHIAW, CHARN, PHINL, PHIBBL, TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUOX, TAUOY
        REAL(rkind)                :: FACTOR, FACTOR2, MWXFINISH, MWYFINISH, EFINISH, DIFF, CONST2, TEMP2
        REAL(rkind)                :: A1BAND, PHIOC, HSTOT, PIBBL, FAGE, FHIGH, FACHFA


!        IF (MINVAL(WALOC) .LT. 0) STOP 'POST 1' 
!        IF (MINVAL(WALOCOLD) .LT. 0) STOP 'POST 2'

        DO IS = 1, NUMSIG
          DO ID = 1, NUMDIR
            AWW3(ID + (IS-1) * NUMDIR)    = WALOC(IS,ID) * CG(IS,IP)
            AWW3OLD(ID + (IS-1) * NUMDIR) = WALOCOLD(IS,ID) * CG(IS,IP) 
          END DO
        END DO

        DO IK=1, NK
          WN2(1+(IK-1)*NTH) = WK(IK,IP)
        END DO
        DO IK=1, NK
          IS0    = (IK-1)*NTH
          DO ITH=2, NTH
            WN2(ITH+IS0) = WN2(1+IS0)
          END DO
        END DO
!
! wind input
!
        AS            = ZERO
        NUMSIG_HF(IP) = NUMSIG
        CALL SET_WIND( IP, WIND10, WINDTH )
        CALL SET_FRICTION( IP, WALOC, WIND10, WINDTH, FPM )
        CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), EMEAN(IP), FMEAN(IP), FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS(IP))

        IF (EMEAN(IP) .LT. THR .AND. WIND10 .GT. THR) CALL SIN_LIN_CAV(IP,WIND10, WINDTH,FPM,SSINL)

        CALL W3SIN4 ( IP, AWW3, CG(:,IP), WN2,  WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, SSINE_WW3, DSSINE_WW3, LLWS, BRLAMBDA)
        CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), EMEAN(IP), FMEAN(IP), FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS(IP))
!
! dissipation 
!
       CALL W3SDS4(AWW3,WK(:,IP),CG(:,IP),UFRIC(IP),USTDIR(IP),DEP(IP),SSDS_WW3,DSSDS_WW3,BRLAMBDA,WHITECAP)
!
! snl4
!
       IF (MESNL .GT. 0) CALL DIASNL4WW3(IP, WNMEAN, WALOC, SSNL4, DSSNL4)
!
! and last but not least the bottom friction dissipation ...
!
       IF (MESBF .GT. 0) CALL SDS_BOTF(IP,WALOC,SSBF,DSSBF)
!
! Resio & Perrie ... scaling ...
!
        FAGE   = FFXFA*TANH(0.3*WIND10*FMEANWS(IP)*PI2/G9)
        FH1    = (FFXFM+FAGE) * FMEAN1
        FH2    = FFXPM / UFRIC(IP) 
        FHIGH  = MIN ( SIG(NK) , MAX ( FH1 , FH2 ) )
        NKH    = MIN ( NK , INT(FACTI2+FACTI1*LOG(MAX(1.E-7,FHIGH))) )
        NKH1   = MIN ( NK , NKH+1 )
        NKH    = MAX ( 2 , MIN ( NKH1 , INT ( FACTI2 + FACTI1*LOG(MAX(1.E-7_rkind,FHIGH)) ) ) )
!
! Add tail
!
        FACHF  = 5.
        FACHFA = XFR**(-FACHF-2)
        DO IK=NKH+1,NK
          DO ITH=1, NTH
            AWW3(ITH+(IK-1)*NTH) = AWW3(ITH+(IK-2)*NTH) * FACHFA + 0. ! Adding a magic zero without a comment ... 
          END DO
        END DO
!
! convert to wwm convetion ...
!
        DO IS = 1, NUMSIG
          DO ID = 1, NUMDIR
            WALOC(IS,ID) = AWW3(ID + (IS-1) * NUMDIR) / CG(IS,IP)
          END DO
        END DO

        CALL CONVERT_VS_VD_WWM(IP, SSINE_WW3, DSSINE_WW3, SSINE, DSSINE)
        CALL CONVERT_VS_VD_WWM(IP, SSDS_WW3, DSSDS_WW3, SSDS, DSSDS)
        IF (MESNL .GT. 0)  CALL CONVERT_WWM_VS_VD(IP, SSNL4_WW3, DSSNL4_WW3, SSNL4, DSSNL4)
        IF (MESBF .GT. 0)  CALL CONVERT_WWM_VS_VD(IP, SSBF_WW3, DSSBF_WW3, SSBF, DSSBF)
!
! compute momentum to BBL 
!
        DO IK=1, NK
          CONST2=DDEN(IK)/CG(IK,IP)*G9/(SIG(IK)/WK(IK,IP))
          DO ITH=1,NTH
            IS=ITH+(IK-1)*NTH
            TEMP2=CONST2*DSSBF_WW3(IS)*AWW3(IS)
            TAUBBL(1) = TAUBBL(1) - TEMP2*ECOS(IS)
            TAUBBL(2) = TAUBBL(2) - TEMP2*ESIN(IS)
          END DO
        END DO
!
! adding the fluxes from waves to ocean ...
!
        WHITECAP(3) = ZERO
        PHIAW       = ZERO
        CHARN       = ZERO
        PHINL       = ZERO
        PHIBBL      = ZERO
        TAUWIX      = ZERO
        TAUWIY      = ZERO
        TAUWNX      = ZERO
        TAUWNY      = ZERO
        HSTOT       = ZERO
       
        DO IK = 1, NK
          FACTOR  = DDEN(IK)/CG(IK,IP)    
          FACTOR2 = FACTOR*G9*WK(IK,IP)/SIG(IK)  
          DO ITH = 1, NTH
            IS      = (IK-1) * NTH + ITH
            PHIAW   = PHIAW  + SSINE_WW3(IS) * DT4S * FACTOR / MAX ( ONE , (ONE-DT4S*DSSINE_WW3(IS))) 
            PIBBL   = PHIBBL - SSBF_WW3(IS) * DT4S  * FACTOR / MAX ( ONE , (ONE-DT4S*DSSBF_WW3(IS))) 
            PHINL   = PHINL  + SSNL4_WW3(IS) * DT4S * FACTOR / MAX ( ONE , (ONE-DT4S*DSSNL4_WW3(IS))) 
            IF (SSINE_WW3(IS) .GT. ZERO) THEN
              WHITECAP(3) = WHITECAP(3) + AWW3(IS) * FACTOR
            ELSE
              ! computes the upward energy flux (counted > 0 upward)
              CHARN = CHARN - (SSINE_WW3(IS))* DT4S * FACTOR / MAX ( ONE , (ONE-DT4S*DSSINE_WW3(IS))) 
            END IF
            HSTOT = HSTOT + AWW3(IS) * FACTOR
          END DO
        END DO
       
        WHITECAP(3) = 4.*SQRT(WHITECAP(3))
        HSTOT  = 4.*SQRT(HSTOT)
        TAUWIX = TAUWIX + TAUWX(IP) * RHOAW * DT4S
        TAUWIY = TAUWIY + TAUWY(IP) * RHOAW * DT4S
        TAUWNX = TAUWNX + TAUWAX * RHOAW * DT4S
        TAUWNY = TAUWNY + TAUWAY * RHOAW * DT4S
!      
!      The wave to ocean flux is the difference between initial energy 
!      and final energy, plus wind input plus the SNL flux to high freq., 
!      minus the energy lost to the bottom boundary layer (BBL) 
!      
        EFINISH  = 0.
        MWXFINISH  = 0.
        MWYFINISH  = 0.
        DO IK = 1, NK
          EBAND = 0.
          A1BAND = 0.
          B1BAND = 0.
          DO ITH = 1, NTH
            DIFF   = AWW3OLD(ITH+(IK-1)*NTH)-AWW3(ITH+(IK-1)*NTH) ! Check that the difference is taken in the right direction!
            EBAND  = EBAND + DIFF
            A1BAND = A1BAND + DIFF*ECOS(ITH)
            B1BAND = B1BAND + DIFF*ESIN(ITH)
          END DO
          EFINISH   = EFINISH  + EBAND * DDEN(IK) / CG(IK,IP)
          MWXFINISH = MWXFINISH  + A1BAND * DDEN(IK) / CG(IK,IP) * WK(IK,IP)/SIG(IK)
          MWYFINISH = MWYFINISH  + B1BAND * DDEN(IK) / CG(IK,IP) * WK(IK,IP)/SIG(IK)
        END DO
!       
!       Transfoe rmation in momentum flux in m^2 / s^2 
!       
        TAUOX = (G9*MWXFINISH+TAUWIX-TAUBBL(1))/DT4S
        TAUOY = (G9*MWYFINISH+TAUWIY-TAUBBL(2))/DT4S
        TAUWIX = TAUWIX/DT4S
        TAUWIY = TAUWIY/DT4S
        TAUWNX = TAUWNX/DT4S
        TAUWNY = TAUWNY/DT4S
!       
!       Transformation in wave energy flux in W/m^2=kg / s^3 
!       
        PHIOC = RHOW*G9*(EFINISH+PHIAW-PHIBBL)/DT4S
        PHIAW = RHOW*G9*PHIAW /DT4S
        PHINL = RHOW*G9*PHINL /DT4S
        PHIBBL= RHOW*G9*PHIBBL/DT4S

#ifdef SCHISM
        OUTT_INTPAR(IP,32) = TAUOX  ! x-component of the wave-ocean momentum flux (tauox in m2.s-2)
        OUTT_INTPAR(IP,33) = TAUOX  ! y-component of the wave-ocean momentum flux (tauoy in m2.s-2) 
        OUTT_INTPAR(IP,34) = PHIOC  ! Wave-to-ocean TKE flux (phioc in W.m-2)
        OUTT_INTPAR(IP,35) = PHIBBL ! phibbl / wave_phibbl : Energy flux due to bottom friction (phioc in W.m-2)
#endif

!        IF (MINVAL(WALOC) .LT. 0) STOP 'POST 3'

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
