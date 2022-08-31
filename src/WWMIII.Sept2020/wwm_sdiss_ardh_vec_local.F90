      SUBROUTINE SDISS_ARDH_VEC_LOCAL (IPP, F, FL, SL, F1MEAN, XKMEAN, &
     &                    PHIEPS, TAUWD, M, SSDS, DSSDS)
! ----------------------------------------------------------------------

!**** *SDISSIP_ARDH_VEC* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.

!     LOTFI AOUF       METEO FRANCE 2013
!     FABRICE ARDHUIN  IFREMER  2013


!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDISSIP (F, FL, IJS, IJL, SL, F1MEAN, XKMEAN,*
!                        PHIEPS, TAUWD, M LCFLX)
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY
!          *F1MEAN* - MEAN FREQUENCY BASED ON 1st MOMENT.
!          *XKMEAN* - MEAN WAVE NUMBER BASED ON 1st MOMENT.
!          *PHIEPS* - ENERGY FLUX FROM WAVES TO OCEAN INTEGRATED OVER 
!                     THE PROGNOSTIC RANGE.
!          *TAUWD*  - DISSIPATION STRESS INTEGRATED OVER
!                     THE PROGNOSTIC RANGE.
!          *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *LCFLX*  - TRUE IF THE CALCULATION FOR THE FLUXES ARE 
!                     PERFORMED.


!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       ARDHUIN et AL. JPO DOI:10.1175/20110JPO4324.1


! ----------------------------------------------------------------------

      !USE YOWFRED  , ONLY : FR, TH, FRATIO, DELTH, GOM, COSTH, SINTH, DFIM
      !USE YOWPCONS , ONLY : RAD     ,G        ,ZPI      ,ROWATER  ,YEPS
      !USE YOWMEAN  , ONLY : EMEAN, FMEAN
      !USE YOWMPP   , ONLY : NINF     ,NSUP
      !USE YOWPARAM , ONLY : NANG     ,NFRE
      !USE YOWSHAL  , ONLY : TFAK     ,INDEP, TCGOND
      !USE YOWSPEC  , ONLY : U10NEW, THWNEW, USNEW
      !USE YOWSTAT  , ONLY : ISHALLO
      !USE YOWTEST  , ONLY : IU06     ,ITEST
      !USE PARKIND1  ,ONLY : JPIM     ,JPRB
      !USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK

      USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, RKIND, &
     &                DFIM, DFIMOFR, DFFR, DFFR2, WK, ZALP, IAB, SWELLFT, &
     &                IUSTAR, IALPHA, USTARM, TAUT, ONETHIRD, RKIND, ONE, &
     &                DELUST, DELALP, BETAMAX, XKAPPA, IDAMPING, TAUWSHELTER, &
     &                FRATIO, EMEAN, USNEW, THWNEW, DEGRAD, LCFLX, MSC, MDC, &
     &                SINTH, COSTH, &
     &                ROWATER => RHOW, &
     &                ROAIR => RHOA, &
     &                TH => SPDIR, &
     &                DELTH => DDIR, &
     &                G => G9, &
     &                ZPI => PI2, &
     &                EPSMIN => SMALL, &
     &                NANG => MDC, &
     &                NFRE => MSC, &
     &                INDEP => DEP, &
     &                CG0 => CG, &
     &                WK0 => WK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: K , M, IPP
      INTEGER :: I, J, I1, J1
      INTEGER :: IS, SDSNTH, DIKCUMUL
      INTEGER :: I_INT, J_INT, M2, K2
      INTEGER :: NSMOOTH(NFRE)
      INTEGER :: ISDSDTH
      INTEGER :: ISSDSBRFDF
      INTEGER :: MIJ
      INTEGER, ALLOCATABLE :: SATINDICES(:,:)

      REAL(rkind) :: XK(NFRE), CG(NFRE)
      REAL(rkind) :: ALFAMEAN
      REAL(rkind) :: TPIINV, TMP00, TMP01, TMP02, TMP03, TMP04  
      REAL(rkind) :: COSWIND
      REAL(rkind) :: DTURB
      REAL(rkind) :: EPSR
      REAL(rkind) :: W, P0, MICHE, DELTA1, DELTA2
      REAL(rkind) :: SCDFM
      REAL(rkind) :: SDISS, YEPS
      REAL(rkind) :: SDSBR, SDSBR2
      REAL(rkind) :: SATURATION2,FACSAT
      REAL(rkind) :: SSDSC3, SSDSC4, SSDSC5, SSDSC6, SDSCOS
      REAL(rkind) :: SSDSHF, SSDSLF, X, DTEMP, TEMP
      REAL(rkind) :: SSDSBRF1, XFR, SXFR,SSDSBR,SDSC1,SSDSP,SSDSC2
      REAL(rkind) :: DSIP_
      REAL(rkind) :: ROG
      REAL(rkind) :: SSDSC1, SSDSC , SSDSISO 
      !REAL(KIND=JPRB) :: ZHOOK_HANDLE
      REAL(rkind) :: SIG(NFRE)
      REAL(rkind) :: BSIGBAJ
      REAL(rkind) :: F1MEAN, XKMEAN, WNMEAN2
      REAL(rkind) :: PHIEPS, TAUWD, CM
      REAL(rkind) :: FACTOR, FACTURB
      REAL(rkind) :: RENEWALFREQ(NANG)
      REAL(rkind), DIMENSION(NFRE) :: CONSTFM
      REAL(rkind) :: BTH0(NFRE)  !saturation spectrum 
      REAL(rkind) :: BTH0S(NFRE)    !smoothed saturation spectrum 
      REAL(rkind) :: BTH(NANG,NFRE)  !saturation spectrum 
      REAL(rkind) :: BTHS(NANG,NFRE)  !smoothed saturation spectrum 
      REAL(rkind),DIMENSION(NANG,NFRE) :: F,FL,SL,A, D
      REAL(rkind),DIMENSION(NANG,NFRE) :: SSDS, DSSDS
      REAL(rkind), ALLOCATABLE, DIMENSION(:) :: C_,C_C,C2_,C2_C2,DSIP_05_C2
      REAL(rkind)   , ALLOCATABLE :: SATWEIGHTS(:,:)
      REAL(rkind)   , ALLOCATABLE :: CUMULW(:,:,:,:)
      LOGICAL :: LLTEST 

      !IF (LHOOK) CALL DR_HOOK('SDISSIP_ARDH_VEC',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

      W=26.
      TPIINV = 1./ZPI
      MIJ = MSC

      ROG = ROWATER*G

      IF (LCFLX) THEN
        PHIEPS = 0.
        TAUWD = 0.
!       !!!! CONSTFM is only defined up to M=MIJ
        SCDFM = 0.5*DELTH*(1.-1./FRATIO)
        DO M=1,MIJ-1
          CONSTFM(M) = ROG*DFIM(M)
        ENDDO
        CONSTFM = ROG*SCDFM*FR(MIJ)
      ENDIF

      XFR = 1.1 !! ??? JEAN BIDLOT: what is XFR is it FRATIO ????
      SSDSBRF1   = 0.5
      SXFR = 0.5*(FRATIO-1/FRATIO)
!      SDSBR     = 1.40E-3       
!      SDSBR     = 1.20E-3       ! from Babanin (personal communication)
      SDSBR     = 9.0E-4
      SDSBR2 = 0.8
      SSDSBR = SDSBR

!      SDSC1  = -2.1 !! This is Bidlot et al. 2005,  Otherwise WAM4 uses -4.5
      SDSC1  = 0.
      SSDSC1  = SDSC1
      SSDSC5  = 0.
      SSDSISO  = 1.
      SSDSP = 2.
!      SSDSC2 = -3.0E-5 ! Retuned by Fabrice from VdW
      SSDSC2 = -2.8E-5
      SSDSC3 = 0.8
      SSDSC4 = 1.
      SSDSC6 = 0.3
      ISDSDTH = 80
      SDSCOS = 2.
      SSDSHF = 1.
      SSDSLF = 1.
      EPSR=SQRT(SSDSBR)

      TMP00 = SSDSC1*ZPI
      YEPS = ROAIR/ROWATER
      TMP01 = SSDSC5/G*YEPS
      DO M=1, NFRE
        SIG(M) = ZPI*FR(M)
      END DO
      ALFAMEAN = (XKMEAN**2)*EMEAN(IPP)
      FACTOR = TMP00     *F1MEAN*(ALFAMEAN*ALFAMEAN)
      FACTURB = TMP01*USNEW(IPP)*USNEW(IPP)
      WNMEAN2 = MAX( 1.E-10 , XKMEAN  )
      IF (ISHALLO.EQ.0) THEN
        DO M=1, NFRE
          XK(M) = WK0(M,IPP)!TFAK(INDEP,M)
          CG(M) = CG0(M,IPP)!TCGOND(INDEP,M)
        END DO
      ELSE
        DO M=1, NFRE
          XK(M) = (SIG(M)**2)/G
          CG(M) = G/(2*ZPI*FR(M))
        END DO
      ENDIF

      DO M=1, NFRE
!cdir outerunroll=4
        DO K=1, NANG
          A(K,M) = TPIINV*CG(M)*F(K,M)
        ENDDO
      END DO


      IF (ISDSDTH.NE.180) THEN
        SDSNTH  = MIN(NINT(ISDSDTH*DEGRAD/(DELTH)),NANG/2-1)
        ALLOCATE(SATINDICES(NANG,SDSNTH*2+1))
        ALLOCATE(SATWEIGHTS(NANG,SDSNTH*2+1))
        DO K=1,NANG
          DO I_INT=K-SDSNTH, K+SDSNTH  
            J_INT=I_INT
            IF (I_INT.LT.1)  J_INT=I_INT+NANG
            IF (I_INT.GT.NANG) J_INT=I_INT-NANG
            SATINDICES(K,I_INT-(K-SDSNTH)+1)=J_INT
            SATWEIGHTS(K,I_INT-(K-SDSNTH)+1)=COS(TH(K)-TH(J_INT))**SDSCOS
          END DO
        END DO
      END IF

!     calcul de cumulw
      SSDSBRF1   = 0.5
      SXFR = 0.5*(FRATIO-1/FRATIO)
! initialise CUMULW
      ALLOCATE(CUMULW(NFRE,NANG,NFRE,NANG))
      DO I=1,NFRE
       DO J=1,NANG
         DO I1=1,NFRE
           DO J1=1,NANG
        CUMULW(I,J,I1,J1)=0.
           ENDDO
         ENDDO
       ENDDO
      END DO
      IF (SSDSC3.NE.0.) THEN

!        DIKCUMUL is the integer difference in frequency bands
!        between the "large breakers" and short "wiped-out waves"
        DIKCUMUL = NINT(SSDSBRF1/(XFR-1.))
        ALLOCATE(C_(NFRE))
        ALLOCATE(C_C(NFRE))
        ALLOCATE(C2_(NFRE-DIKCUMUL))
        ALLOCATE(C2_C2(NFRE-DIKCUMUL))
        ALLOCATE(DSIP_05_C2(NFRE-DIKCUMUL))
        DO M=1,NFRE  
          C_(M)=G/SIG(M)  ! Valid in deep water only
          C_C(M)=C_(M)*C_(M)
        END DO
        DO M=1,NFRE  
          DO K=1,NANG
        DO M2=1,M-DIKCUMUL
          C2_(M2)=G/SIG(M2)
          C2_C2(M2)=C2_(M2)*C2_(M2)
          DSIP_ = SIG(M2)*SXFR
          DSIP_05_C2(M2)=DSIP_/(0.5*C2_(M2))
            DO K2=1,NANG
                CUMULW(M,K,M2,K2)=SQRT(C_C(M)+C2_C2(M2)    &
     &          -2*C_(M)*C2_(M2)*COSTH(1+ABS(K2-K)))*DSIP_05_C2(M2) 
              END DO
            END DO 
          END DO
        END DO
        DEALLOCATE(C_)
        DEALLOCATE(C_C)
        DEALLOCATE(C2_)
        DEALLOCATE(C2_C2)
        DEALLOCATE(DSIP_05_C2)

! Multiplies by lambda(k,theta)=1/(2*pi**2) and 
! and the coefficient that transforms  SQRT(B) to Banner et al. (2000)'s epsilon
! 2.26 is equal to 5.55 (Banner & al. 2000) times 1.6**2 / 2pi where
! 1.6 is the ratio between Banner's epsilon and SQRT(B)

        TMP02 = 2*TPIINV*2.26
        DO I=1,NFRE
         DO J=1,NANG 
           DO I1=1,NFRE
             DO J1=1,NANG
               CUMULW(I,J,I1,J1)=CUMULW(I,J,I1,J1)*TMP02
             END DO
           ENDDO
         END DO
        END DO
      END IF
        DO  M=1, NFRE
          FACSAT=(XK(M)**3)*DELTH
          BTH0(M)=SUM(A(1:NANG,M))*FACSAT
        END DO
!      DO K=1,NANG
      DO  M=1, NFRE
        DO K=1,NANG
          FACSAT=(XK(M)**3)*DELTH
          ! integrates around full circle
          BTH(K,M) = SUM(SATWEIGHTS(K,1:SDSNTH*2+1)*A(SATINDICES(K,1:SDSNTH*2+1),M))*FACSAT
        END DO
        BTH0(M) = MAXVAL(BTH(1:NANG,M))
      END DO

!/ST3      SDSBR     = 1.20E-3 ! Babanin (personnal communication)
      ISSDSBRFDF  = 22    ! test pour DC
      ISSDSBRFDF  = 0
!/ST3      SDSBRF1   = 0.5
!/ST3      SDSBRF2   = 0.
      IF (ISSDSBRFDF.GT.0.AND.ISSDSBRFDF.LT.NFRE/2) THEN 
!cdir collapse
        BTH0S=BTH0
!cdir collapse
        NSMOOTH=ONE
!cdir collapse
        BTHS=BTH
!cdir outerunroll=4
        DO M=1, ISSDSBRFDF
          BTH0S  (1+ISSDSBRFDF)=BTH0S  (1+ISSDSBRFDF)+BTH0(M)
          NSMOOTH(1+ISSDSBRFDF)=NSMOOTH(1+ISSDSBRFDF)+1
!cdir collapse
          DO K=1,NANG       
            BTHS(K,M)=BTHS(K,M)+BTH(K,M)
          END DO 
        ENDDO
        DO M=2+ISSDSBRFDF,1+2*ISSDSBRFDF
!cdir nodep
          BTH0S  (1+ISSDSBRFDF)=BTH0S  (1+ISSDSBRFDF)+BTH0(M)
          NSMOOTH(1+ISSDSBRFDF)=NSMOOTH(1+ISSDSBRFDF)+1
!cdir collapse
          DO K=1,NANG       
            BTHS(K,M)=BTHS(K,M)+BTH(K,M)
          END DO
        ENDDO
        DO M=ISSDSBRFDF,1,-1
!cdir nodep
          BTH0S  (M)=BTH0S  (M+1)-BTH0(M+ISSDSBRFDF+1)
          NSMOOTH(M)=NSMOOTH(M+1)-1
!cdir collapse
          DO K=1,NANG
            BTHS(K,M)=BTHS(K,M)-BTH(K,M)
          END DO
        ENDDO
        DO M=2+ISSDSBRFDF,NFRE-ISSDSBRFDF
!cdir nodep
          BTH0S  (M)=BTH0S  (M-1)-BTH0(M-ISSDSBRFDF-1)+BTH0(M+ISSDSBRFDF)
          NSMOOTH(M)=NSMOOTH(M-1)
!cdir collapse
          DO K=1,NANG       
            BTHS(K,M)=BTHS(K,M)-BTH(K,M)+BTH(K,M)
          END DO
        ENDDO
!cdir novector
        DO M=NFRE-ISSDSBRFDF+1,NFRE
!cdir nodep
          BTH0S  (M)=BTH0S  (M-1)-BTH0(M-ISSDSBRFDF)
          NSMOOTH(M)=NSMOOTH(M-1)-1
!cdir collapse
          DO K=1,NANG       
            BTHS(K,M)=BTHS(K,M)-BTH(K,M)
          END DO
        END DO

!  final division by NSMOOTH

!cdir collapse
        DO M=1,NFRE
          BTH0(M)=MAX(0.,BTH0S(M)/NSMOOTH(M))
        END DO 
        DO M=1,NFRE
!cdir outerunroll=4
          DO K=1,NANG
            BTH(K,M)=MAX(0.,BTHS(K,M)/NSMOOTH(M))
          END DO
        END DO 
           
      END IF

!      DELTA1 = 0.4
!      DELTA2 = 0.6
      DELTA1 = 0.
      DELTA2 = 0.
      MICHE = 1.0
      TMP03 = 1.0/(SSDSBR*MICHE)

      DO  M=1, NFRE

        LLTEST = (SSDSC3.NE.0.AND.M.GT.DIKCUMUL)

        IF (XKMEAN.NE.0) THEN
          X           = WK0(M,IPP)/XKMEAN!TFAK(INDEP(IJ),M)/XKMEAN(IJ)
          BSIGBAJ = FACTOR*( (1.-DELTA2)*X + DELTA2*X**2)
        ELSE
          BSIGBAJ = 0
        ENDIF

        IF (ISHALLO.EQ.0) THEN
          CM=WK0(M,IPP)/SIG(M)!TFAK(INDEP(IJ),M)/SIG(M)
        ELSE
!AR: below is a nice bug always the last index ... in the original, should be IJ
          CM=SIG(M)/G
        ENDIF

        DO K=1,NANG

          RENEWALFREQ(K)=0.
          ! Correction of saturation level for shallow-water kinematics
          ! Cumulative effect based on lambda   (breaking probability is
          ! the expected rate of sweeping by larger breaking waves)
          IF (LLTEST) THEN
            DO M2=1,M-DIKCUMUL  
              !AR: below the index is wrong ...
                DO K2=1,NANG
                    IF (BTH0(M2).GT.SSDSBR) THEN
                  ! Integrates over frequencies M2 and directions K2 to 
                  ! Integration is performed from M2=1 to a frequency lower than M: IK-DIKCUMUL
                 RENEWALFREQ(K)=RENEWALFREQ(K)+ CUMULW(M,K,M2,K2) &
     &              *(MAX(SQRT(BTH(K2,M2))-EPSR,0.))**2
                    ENDIF
                END DO
            END DO
          ENDIF

          SATURATION2=TANH(10*(((BTH(K,M)/SSDSBR)**0.5)-SDSBR2))
          COSWIND=(COSTH(K)*COS(THWNEW(IPP))+SINTH(K)*SIN(THWNEW(IPP)))   ! ÃvÃ©rifier K ?
          DTURB=-2.*SIG(M)*XK(M)*FACTURB*COSWIND  ! Theory -> stress direction
          P0=SSDSP ! -0.5*SSDSC3*(1-TANH(W*USTAR*XK(M)/SIG(M)-0.1))  ! for SDSC3=1 this is vdW et al. 

          TMP04 = SSDSC3*RENEWALFREQ(K)
!          DTEMP=SSDSC2 * SIG(M) &
          DTEMP=SSDSC2 * SIG(M) &
     &  * (  SSDSC6 *(MAX(0.,BTH0(    M)*TMP03-SSDSC4))**P0 &
     &  + (1-SSDSC6)*(MAX(0.,BTH (K,M)*TMP03-SSDSC4))**P0)&
     &  - (TMP04+DTURB)  !terme cumulatif
          D(K,M) = DTEMP + BSIGBAJ*SSDSLF *0.5*(1-SATURATION2) &
     &                  + BSIGBAJ*SSDSHF *0.5*(SATURATION2+1)
          WRITE(111116,'(10F15.8)') D(K,M),DTEMP,SIG(M),RENEWALFREQ(K),BTH0(M),BTH(K,M)
        END DO

!cdir outerunroll=4
        DO K=1, NANG
          SL(K,M) = SL(K,M)+D(K,M)*F(K,M)
          FL(K,M) = FL(K,M)+D(K,M)
          SSDS(K,M) = D(K,M)*F(K,M)
          DSSDS(K,M) = D(K,M)
          IF (LCFLX.AND.M.LE.MIJ) THEN
            SDISS = D(K,M)*F(K,M)
            PHIEPS = PHIEPS+SDISS*CONSTFM(M)
            TAUWD  = TAUWD+CM*SDISS*CONSTFM(M)
          ENDIF
        END DO
      END DO

      IF (ALLOCATED(CUMULW)) DEALLOCATE (CUMULW)
      IF (ALLOCATED(SATWEIGHTS)) DEALLOCATE (SATWEIGHTS)
      IF (ALLOCATED(SATINDICES)) DEALLOCATE (SATINDICES)

      !IF (LHOOK) CALL DR_HOOK('SDISSIP_ARDH_VEC',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE SDISS_ARDH_VEC_LOCAL
