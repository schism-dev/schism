#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SOURCE_INT_EXP
         USE DATAPOOL
#ifdef ST41
         USE W3SRC4MD_OLD, ONLY : LFIRSTSOURCE
#elif ST42
         USE W3SRC4MD, ONLY : LFIRSTSOURCE
#endif
         IMPLICIT NONE

         INTEGER        :: IP, IS, ID, I, K, M
         INTEGER        :: NIT_SIN, NIT_SDS, NIT_SNL4, NIT_SNL3, NIT_SBR, NIT_SBF, NIT_ALL
         INTEGER, SAVE  :: IFIRST, ISELECT
         REAL(rkind)    :: ACLOC(MSC,MDC), IMATRA(MSC,MDC), IMATDA(MSC,MDC), SSBRL2(MSC,MDC)
         REAL(rkind)    :: DT4S_T, DT4S_E, DT4S_Q, DT4S_H, DT4S_TQ, DT4S_TS, VEC2RAD

!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP&         SHARED(MNP,MSC,MDC,DEP,DMIN,IOBP,SMETHOD, &
!$OMP&         DT4S,DTMIN_SDS,DTMIN_SIN,DTMIN_SBR,DTMIN_DYN, &
!$OMP&         DTMIN_SNL3, DTMIN_SNL4, DTMIN_SBF, WINDXY, WINDFAC,&
!$OMP&         LSOUBOUND,ISHALLOW,LADVTEST,LMAXETOT, &
!$OMP&         NDYNITER_SIN,NDYNITER_SNL4, NDYNITER_SDS, &
!$OMP&         NDYNITER_SBR, NDYNITER_SNL3, NDYNITER_SBF, &
!$OMP&         NDYNITER,LSOURCESWAM,LLIMT,MESTR,MESBR,MESBF,AC2,AC1,&
!$OMP&         FL,SL,FL3,THWNEW,THWOLD,U10NEW,U10OLD,Z0NEW,Z0OLD) &
!$OMP&         PRIVATE(IP,IS,ID,ACLOC,NIT_SIN,NIT_SDS,&
!$OMP&         NIT_SNL4,NIT_SNL3,NIT_SBR,NIT_SBF,NIT_ALL,SSBRL2,&
!$OMP&         IMATRA,IMATDA,VEC2RAD)
!$OMP DO SCHEDULE(DYNAMIC,1)
         DO IP = 1, MNP
!           IF (IP_IS_STEADY(IP) .EQ. 1) CYCLE
           IF (DEP(IP) .LT. DMIN) CYCLE
           IF (IOBP(IP) .EQ. 0) THEN
               ACLOC  = AC2(:,:,IP)
               IF (SMETHOD == 1) THEN
                 IF (MESBR .GT. 0 .OR. MESBF .GT. 0) CALL RKS_SP3(IP,DT4S,.FALSE.,ACLOC,30)
                 IF (LSOURCESWAM) THEN
                   CALL SOURCE_INT_EXP_WAM(IP, ACLOC)  
                 ELSE
                   CALL INT_IP_STAT(IP,DT4S,LLIMT,ACLOC,20)
                 ENDIF
                 IF (MESTR .GT. 0) CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_DYN, NDYNITER, ACLOC, NIT_ALL, 4)
               ELSE IF (SMETHOD == 2) THEN
                 CALL INT_IP_STAT(IP,DT4S, LLIMT,ACLOC, 10)
               ELSE IF (SMETHOD == 3) THEN
                 CALL RKS_SP3(IP,DT4S,LLIMT,ACLOC, 10)
               ELSE IF (SMETHOD == 4) THEN
                 CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_DYN, NDYNITER, ACLOC, NIT_ALL, 10)
               ELSE IF (SMETHOD == 5) THEN ! Full splitting of all source embedded within a dynamic RK-3 Integration ... 
                 CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SIN,  NDYNITER_SIN  , ACLOC, NIT_SIN, 1) ! Sin
                 CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SNL4,  NDYNITER_SNL4 , ACLOC, NIT_SNL4, 2)! Snl4b
                 CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SDS,  NDYNITER_SDS  , ACLOC, NIT_SDS, 3) ! Sds
                 CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SNL3,  NDYNITER_SNL3 , ACLOC, NIT_SNL3, 4)! Snl3
                 CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SBR,  NDYNITER_SBR  , ACLOC, NIT_SBR, 5) ! Sbr
                 CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SBF,  NDYNITER_SBF  , ACLOC, NIT_SBF, 6) ! Sbf
               END IF
               CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .TRUE.,1,'RECALC DOMAIN') ! Update everything based on the new spectrum ...
               IF (LMAXETOT .AND. .NOT. LADVTEST .AND. ISHALLOW(IP) .EQ. 1) THEN
                 CALL BREAK_LIMIT(IP,ACLOC,SSBRL2)
               ENDIF
               AC2(:,:,IP) = ACLOC
           ELSE !Boundary node ... 
             IF (LSOUBOUND) THEN ! Source terms on boundary ...
               IF (IOBP(IP) .NE. 2) THEN
                 ACLOC  = AC2(:,:,IP)
                 IF (SMETHOD == 1) THEN
                   CALL RKS_SP3(IP,DT4S,.FALSE.,ACLOC,30)
                   IF (LSOURCESWAM) THEN
                     CALL SOURCE_INT_EXP_WAM(IP, ACLOC)
                   ELSE
                     CALL INT_IP_STAT(IP,DT4S,LLIMT,ACLOC,20)
                   ENDIF
                   CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_DYN, NDYNITER, ACLOC, NIT_ALL,4)
                 ELSE IF (SMETHOD == 2) THEN
                   CALL INT_IP_STAT(IP,DT4S, LLIMT,ACLOC,10)
                 ELSE IF (SMETHOD == 3) THEN
                   CALL RKS_SP3(IP,DT4S,LLIMT,ACLOC,10)
                 ELSE IF (SMETHOD == 4) THEN
                   CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_DYN, NDYNITER, ACLOC, NIT_ALL,10)
                 ELSE IF (SMETHOD == 5) THEN ! Full splitting of all source embedded within a dynamic RK-3 Integration ... 
                   CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SIN,  NDYNITER_SIN  , ACLOC, NIT_SIN,1) ! Sin
                   CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SNL4,  NDYNITER_SNL4 , ACLOC, NIT_SNL4,2)! Snl4b
                   CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SDS,  NDYNITER_SDS  , ACLOC, NIT_SDS,3) ! Sds
                   CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SNL3,  NDYNITER_SNL3 , ACLOC, NIT_SNL3,4)! Snl3
                   CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SBR,  NDYNITER_SBR  , ACLOC, NIT_SBR,5) ! Sbr
                   CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SBF,  NDYNITER_SBF  , ACLOC, NIT_SBF,6) ! Sbf
                 END IF
                 CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .TRUE., 1, 'RECALC BOUNDARY') ! Update everything based on the new spectrum ...
                 IF (LMAXETOT .AND. .NOT. LADVTEST .AND. ISHALLOW(IP) .EQ. 1) THEN
                   CALL BREAK_LIMIT(IP,ACLOC,SSBRL2)
                 ENDIF
                 AC2(:,:,IP) = ACLOC
               ENDIF
             ELSE
               ACLOC = AC2(:,:,IP)
               IF (LMAXETOT .AND. .NOT. LADVTEST .AND. ISHALLOW(IP) .EQ. 1) THEN ! limit wave height on the boundary ...
                 CALL BREAK_LIMIT(IP,ACLOC,SSBRL2)
                 AC2(:,:,IP) = ACLOC
               ENDIF
             ENDIF
           ENDIF
           AC1(:,:,IP) = AC2(:,:,IP)
         ENDDO
!$OMP END PARALLEL 

#if defined ST41 || defined ST42
         LFIRSTSOURCE = .FALSE.
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SOURCE_INT_IMP_WWM
         USE DATAPOOL
#ifdef ST41
         USE W3SRC4MD_OLD, ONLY : LFIRSTSOURCE
#endif
#ifdef ST42
         USE W3SRC4MD, ONLY : LFIRSTSOURCE
#endif
         IMPLICIT NONE

         INTEGER :: IP

         REAL(rkind)    :: ACLOC(MSC,MDC)
         REAL(rkind)    :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

!$OMP WORKSHARE
         IMATDAA = 0.
         IMATRAA = 0.
!$OMP END WORKSHARE

!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP&         SHARED(MNP,MSC,MDC,DEP,DMIN,IOBP,SMETHOD, &
!$OMP&         LSOUBOUND) &
!$OMP&         PRIVATE(IP,ACLOC,AC2,IMATDA,IMATRA,IMATRAA,IMATDAA)
!$OMP DO SCHEDULE(DYNAMIC,1)
         DO IP = 1, MNP
           IF (DEP(IP) .LT. DMIN) CYCLE
           IF ((ABS(IOBP(IP)) .NE. 1 .AND. IOBP(IP) .NE. 3)) THEN
             IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
               ACLOC = AC2(:,:,IP)
               CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE., 10, 'SOURCE_INT_IMP_WWM DOMAIN') 
               !CALL CYCLE3(IP, ACLOC, IMATRA, IMATDA)
               IMATDAA(:,:,IP) = IMATDA 
               IMATRAA(:,:,IP) = IMATRA 
             ENDIF
           ELSE
             IF (LSOUBOUND) THEN ! Source terms on boundary ...
               IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
                 ACLOC = AC2(:,:,IP)
                 CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE.,10, 'SOURCE_INT_IMP_WWM BOUNDARY')
                 !CALL CYCLE3(IP, ACLOC, IMATRA, IMATDA)
                 IMATDAA(:,:,IP) = IMATDA
                 IMATRAA(:,:,IP) = IMATRA 
               ENDIF
             ENDIF
           ENDIF
         END DO
!$OMP END PARALLEL 
#if defined ST41 || defined ST42
         LFIRSTSOURCE = .FALSE.
#endif

!         PAUSE
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************      
      SUBROUTINE SOURCE_INT_EXP_WAM(IP, ACLOC)
         USE DATAPOOL
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: IP
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)
         
         INTEGER      :: IS, ID, IMETHOD, ICODE
         REAL(rkind)  :: VEC2RAD
         REAL(rkind)  :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind)  :: ETOT,SME01,SME10,KME01,KMWAM,KMWAM2,HS,WIND10
         REAL(rkind)  :: ETAIL,EFTAIL,EMAX,LIMAC,NEWDAC,FPM,WINDTH,TEMP,GTEMP1
         REAL(rkind)  :: RATIO,LIMFAC,LIMDAC,GTEMP2,FLHAB,DELFL,USFM, NEWDACDT
         REAL(rkind)  :: MAXDAC, MAXDACDT, MAXDACDTDA, SC, SP, DNEWDACDTDA, JAC, FF

         REAL(rkind),DIMENSION(MDC,MSC)  :: SSDS,DSSDS,SSNL4,DSSNL4,SSIN,DSSIN
     
         ICODE = 1 
         
         IF (MESIN .GT. 0 .OR. MESDS .GT. 0 .OR. MESNL .GT. 0) THEN
           DO IS = 1, MSC
             DO ID = 1, MDC
               FL3(IP,ID,IS) = ACLOC(IS,ID) * PI2 * SPSIG(IS)
               FL(IP,ID,IS)  = FL3(IP,ID,IS)
               SL(IP,ID,IS)  = FL(IP,ID,IS)
             END DO
           END DO
           THWOLD(IP,1) = THWNEW(IP)
           U10NEW(IP) = MAX(TWO,SQRT(WINDXY(IP,1)**2+WINDXY(IP,2)**2))*WINDFAC
           Z0NEW(IP) = Z0OLD(IP,1)
           THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))
           IF (ICODE == 1) THEN
             CALL IMPLSCH (FL3(IP,:,:), FL(IP,:,:), IP, IP, 1, &
     &                     THWOLD(IP,1), USOLD(IP,1), &
     &                     TAUW(IP), Z0OLD(IP,1), &
     &                     ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                     U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                     Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                     SL(IP,:,:), FCONST(IP,:))
           ELSE IF (ICODE == 2) THEN
             CALL PREINTRHS (FL3(IP,:,:), FL(IP,:,:), IP, IP, 1, &
     &                       THWOLD(IP,1), USOLD(IP,1), &
     &                       TAUW(IP), Z0OLD(IP,1), &
     &                       ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                       U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                       Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                       SL(IP,:,:), FCONST(IP,:), FMEANWS(IP), MIJ(IP), &
     &                       SSDS, DSSDS, SSIN, DSSIN, &
     &                       SSNL4, DSSNL4)
             CALL INTSPECWAM (FL3(IP,:,:), FL(IP,:,:), IP, IP, 1, &
     &                      THWOLD(IP,1), USOLD(IP,1), &
     &                      TAUW(IP), Z0OLD(IP,1), &
     &                      ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                      U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                      Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                      SL(IP,:,:), FCONST(IP,:), FMEANWS(IP), MIJ(IP))
             CALL POSTINTRHS (FL3(IP,:,:), FL(IP,:,:), IP, IP, 1, &
     &                      THWOLD(IP,1), USOLD(IP,1), &
     &                      TAUW(IP), Z0OLD(IP,1), &
     &                      ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                      U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                      Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                      SL(IP,:,:), FCONST(IP,:), FMEANWS(IP), MIJ(IP))
           ELSE IF (ICODE == 3) THEN
             CALL IMPLSCH_LOCAL (IP, FL3(IP,:,:), FL(IP,:,:), 1, &
     &                     SL(IP,:,:))
           ENDIF ! ICODE
           DO IS = 1, MSC
             DO ID = 1, MDC
               ACLOC(IS,ID) =  MAX(ZERO, FL3(IP,ID,IS) / PI2 / SPSIG(IS))
             END DO
           END DO
         ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SOURCE_INT_IMP_WAM_PRE
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER      :: IP, IS, ID, IMETHOD
         REAL(rkind)  :: ACLOC(MSC,MDC), VEC2RAD
         REAL(rkind)  :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind)  :: SSNL3(MSC,MDC),DSSNL3(MSC,MDC), DSSVEG(MSC,MDC)
         REAL(rkind)  :: SSBR(MSC,MDC),DSSBR(MSC,MDC), SSVEG(MSC,MDC)
         REAL(rkind)  :: SSBF(MSC,MDC),DSSBF(MSC,MDC), SSINL(MSC,MDC)
         REAL(rkind)  :: ETOT,SME01,SME10,KME01,KMWAM,KMWAM2,HS,WIND10
         REAL(rkind)  :: ETAIL,EFTAIL,EMAX,LIMAC,NEWDAC,FPM,WINDTH,TEMP,GTEMP1
         REAL(rkind)  :: RATIO,LIMFAC,LIMDAC,GTEMP2,FLHAB,DELFL,USFM, NEWDACDT
         REAL(rkind)  :: MAXDAC, MAXDACDT, MAXDACDTDA, SC, SP, DNEWDACDTDA, JAC, FF

         REAL(rkind),DIMENSION(MDC,MSC)  :: SSDS,DSSDS,SSNL4,DSSNL4,SSIN,DSSIN

         IMETHOD = 4 

!$OMP WORKSHARE
         IMATDAA = ZERO
         IMATRAA = ZERO
!$OMP END WORKSHARE

!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP&         SHARED(MNP,MSC,MDC,DEP,DMIN,IOBP,SMETHOD, &
!$OMP&         LSOUBOUND,LINID,MESTR,MESBR,MESBF,LIMFAK,&
!$OMP&         COFRM4,WINDFAC,MESIN,MESDS,MESNL,SPSIG) &
!$OMP&         PRIVATE(IP,ACLOC,AC2,IMATDA,IMATRA,IMATRAA,IMATDAA, &
!$OMP&         SSNL3,DSSNL3,SSBR,DSSBR,SSBF,DSSBF,FL3,FL,Z0NEW,&
!$OMP&         SSDS,DSSDS,JAC,USNEW,ZIDLNEW,FF,FMEAN,USFM,LIMFAC,&
!$OMP&         FCONST,ZIDLOLD,TAUW,USOLD,WINDXY,WIND10,WINDTH,&
!$OMP&         TEMP,HS,ETOT,FPM,ISHALLOW,KME01,SME10,KMWAM, &
!$OMP&         NEWDAC,KMWAM2,MAXDAC,WK,SC,CG,SME01,GTEMP2,&
!$OMP&         DELFL,FLHAB,GTEMP1,IMETHOD,DT4A,SSNL4,DSSNL4,&
!$OMP&         SSIN,DSSIN,FMEANWS,ROAIRN,MIJ,Z0OLD,ROAIRO, &
!$OMP&         THWOLD,THWNEW,U10NEW,SL,SSINL)
!$OMP DO SCHEDULE(DYNAMIC,1)
         DO IP = 1, MNP
           SSNL3 = ZERO; DSSNL3 = ZERO
           SSBR  = ZERO; DSSBR  = ZERO
           SSBF  = ZERO; DSSBF  = ZERO
           SSINL = ZERO
           IF (DEP(IP) .LT. DMIN) THEN
             IMATRAA(:,:,IP) = ZERO
             IMATDAA(:,:,IP) = ZERO
             CYCLE
           ENDIF
           IF (MESIN .GT. 0 .OR. MESDS .GT. 0 .OR. MESNL .GT. 0) THEN
             DO IS = 1, MSC
               DO ID = 1, MDC
                 FL3(IP,ID,IS) = AC2(IS,ID,IP) * PI2 * SPSIG(IS)
                 FL(IP,ID,IS)  = FL3(IP,ID,IS)
                 SL(IP,ID,IS)  = FL(IP,ID,IS)
               END DO
             END DO
             THWOLD(IP,1) = THWNEW(IP)
             U10NEW(IP) = MAX(TWO,SQRT(WINDXY(IP,1)**2+WINDXY(IP,2)**2))*WINDFAC
             Z0NEW(IP) = Z0OLD(IP,1)
             THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))
             CALL PREINTRHS (FL3(IP,:,:), FL(IP,:,:), IP, IP, 1, &
     &                       THWOLD(IP,1), USOLD(IP,1), &
     &                       TAUW(IP), Z0OLD(IP,1), &
     &                       ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                       U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                       Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                       SL(IP,:,:), FCONST(IP,:), FMEANWS(IP), MIJ(IP), &
     &                       SSDS, DSSDS, SSIN, DSSIN, &
     &                       SSNL4, DSSNL4)
           ENDIF ! MESIN .EQ. 0 .AND. MESDS .AND. 0 .AND. MESNL = 0
           IF (IOBP(IP) .EQ. 0) THEN
             DO ID = 1, MDC
               DO IS = 1, MSC 
                 IF (AC2(IS,ID,IP) .LT. THR) CYCLE
                 JAC = ONE/PI2/SPSIG(IS)
                 !IMATDAA(IS,ID,IP) =  FL(IP,ID,IS) !... this is not working right, reason is unknown, there signchanges that should not be
                 !IMATRAA(IS,ID,IP) =  SL(IP,ID,IS)/PI2/SPSIG(IS) 
                 FF = FL3(IP,ID,IS)
                 !WRITE(11140,'(2I10,7E20.10)') IS, ID, FF, SSDS(ID,IS)/FF, DSSDS(ID,IS), SSIN(ID,IS)/FF, DSSNL4(ID,IS), SSNL4(ID,IS)/FF, DSSNL4(ID,IS) 
                 IF (IMETHOD == 0) THEN
                   IMATRAA(IS,ID,IP) =  SL(IP,ID,IS)/PI2/SPSIG(IS)
                   IMATDAA(IS,ID,IP) = ZERO 
                 ELSE IF (IMETHOD == 1) THEN 
                   IMATRAA(IS,ID,IP) = (SSIN(ID,IS)+SSDS(ID,IS)+SSNL4(ID,IS))*JAC
                 ELSE IF (IMETHOD == 2) THEN
                   IMATRAA(IS,ID,IP) = (SSIN(ID,IS)+SSNL4(ID,IS))*JAC
                   IMATDAA(IS,ID,IP) = -TWO*DSSDS(ID,IS)
                 ELSE IF (IMETHOD == 3) THEN
                   IF (SSIN(ID,IS) .GT. ZERO) THEN
                     IMATRAA(IS,ID,IP) = SSIN(ID,IS)*JAC
                   ELSE
                     IMATDAA(IS,ID,IP) = -TWO*DSSIN(ID,IS)
                   ENDIF
                 ELSE IF (IMETHOD == 4) THEN
                   GTEMP1 = MAX((1.-DT4A*FL(IP,ID,IS)),1.)
                   GTEMP2 = DT4A*SL(IP,ID,IS)/GTEMP1
                   FLHAB  = ABS(GTEMP2)
                   DELFL  = COFRM4(IS)*DT4A
                   USFM   = USNEW(IP)*MAX(FMEANWS(IP),FMEAN(IP))
                   TEMP   = USFM*DELFL
                   FLHAB  = MIN(FLHAB,TEMP)
                   IMATRAA(IS,ID,IP) = SIGN(FLHAB,GTEMP2)/DT4A*JAC
                   LIMFAC            = MIN(ONE,ABS(SIGN(FLHAB,GTEMP2/DT4A))/MAX(SMALL,ABS(IMATRAA(IS,ID,IP))))
                   IMATDAA(IS,ID,IP) = ZERO ! -FL(IP,ID,IS)
                 ENDIF 
               ENDDO ! ID
             ENDDO ! IS
             IF (LNANINFCHK) THEN
               IF (SUM(IMATRAA(:,:,IP)) .NE. SUM(IMATRAA(:,:,IP))) THEN
                 WRITE(12,*) 'NAN AT NODE', IP, 'AT DEPTH', DEP(IP)
                 CALL WWM_ABORT('NAN IN SOURCES 1a') 
               ENDIF
               IF (SUM(IMATDAA(:,:,IP)) .NE. SUM(IMATDAA(:,:,IP))) THEN
                 WRITE(12,*) 'NAN AT NODE', IP, 'AT DEPTH', DEP(IP)
                 CALL WWM_ABORT('NAN IN SOURCES 1b') 
               ENDIF
             ENDIF
             ACLOC = AC2(:,:,IP)
             IF (.NOT. LINID) THEN
               CALL SET_WIND( IP, WIND10, WINDTH )
               CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
               CALL SIN_LIN_CAV(IP,WIND10,WINDTH,FPM,IMATRA,SSINL)
               IMATRAA(:,:,IP) = IMATRAA(:,:,IP) + SSINL
             ENDIF
             IF (ISHALLOW(IP) .EQ. 1) THEN
               CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)
               IF (MESTR .GT. 0) THEN
                 CALL TRIAD_ELDEBERKY(IP, HS, SME01, ACLOC, IMATRA, IMATDA, SSNL3, DSSNL3)
                 DO IS = 1, MSC
                   DO ID = 1, MDC
                     NEWDAC = SSNL3(IS,ID)*DT4A/MAX((1.-DT4A*DSSNL3(IS,ID)),1.)
                     MAXDAC = 0.0081*LIMFAK/(TWO*SPSIG(IS)*WK(IS,IP)**3*CG(IS,IP))*100
                     LIMFAC = ONE/MAX(ONE,NEWDAC/MAXDAC)
                     SC = SIGN(MIN(ABS(NEWDAC),MAXDAC),NEWDAC)/DT4A
                     !SSNL3(IS,ID)  = SC
                     !DSSNL3(IS,ID) = DSSNL3(IS,ID)*LIMFAC
                     !IF (ABS(SC) .GT. THR) WRITE(*,'(2I10,5F20.8)') IS, ID, NEWDAC, MAXDAC, DSSNL3(IS,ID), LIMFAC
                   END DO
                 END DO
               ENDIF ! MESTR
               IF (MESBR .GT. 0) CALL SDS_SWB(IP, SME01, KMWAM, ETOT, HS, ACLOC, IMATRA, IMATDA, SSBR, DSSBR)
               IF (MESBF .GT. 0) CALL SDS_BOTF(IP,ACLOC,IMATRA,IMATDA,SSBF,DSSBF)
               IMATDAA(:,:,IP) = IMATDAA(:,:,IP) + DSSBR  + DSSNL3 + DSSBF + DSSVEG 
               IMATRAA(:,:,IP) = IMATRAA(:,:,IP) + SSBR + SSNL3 + SSVEG 
             ENDIF ! ISHALLOW(IP) .EQ. 1
           ELSE ! IOBP(IP) .NE. 0
             IF (LSOUBOUND) THEN ! Source terms on boundary ...
               IF (IOBP(IP) .NE. 2) THEN
                 DO ID = 1, MDC
                   DO IS = 1, MSC
                     IF (AC2(IS,ID,IP) .LT. THR) CYCLE
                     JAC = ONE/PI2/SPSIG(IS)
                     !IMATDAA(IS,ID,IP) =  FL(IP,ID,IS) !... this is not working right, reason is unknown, there signchanges that should not be
                     !IMATRAA(IS,ID,IP) =  SL(IP,ID,IS)/PI2/SPSIG(IS) 
                     FF = FL3(IP,ID,IS)
                     !WRITE(11140,'(2I10,7E20.10)') IS, ID, FF, SSDS(ID,IS)/FF, DSSDS(ID,IS), SSIN(ID,IS)/FF, DSSNL4(ID,IS), SSNL4(ID,IS)/FF, DSSNL4(ID,IS) 
                     IF (IMETHOD == 0) THEN
                       IMATRAA(IS,ID,IP) =  SL(IP,ID,IS)/PI2/SPSIG(IS)
                       IMATDAA(IS,ID,IP) = ZERO
                     ELSE IF (IMETHOD == 1) THEN
                       IMATRAA(IS,ID,IP) = (SSIN(ID,IS)+SSDS(ID,IS)+SSNL4(ID,IS))*JAC
                     ELSE IF (IMETHOD == 2) THEN
                       IMATRAA(IS,ID,IP) = (SSIN(ID,IS)+SSNL4(ID,IS))*JAC
                       IMATDAA(IS,ID,IP) = -TWO*DSSDS(ID,IS)
                     ELSE IF (IMETHOD == 3) THEN
                       IF (SSIN(ID,IS) .GT. ZERO) THEN
                         IMATRAA(IS,ID,IP) = SSIN(ID,IS)*JAC
                       ELSE
                         IMATDAA(IS,ID,IP) = -TWO*DSSIN(ID,IS)
                       ENDIF
                     ELSE IF (IMETHOD == 4) THEN
                       GTEMP1 = MAX((1.-DT4A*FL(IP,ID,IS)),1.)
                       GTEMP2 = DT4A*SL(IP,ID,IS)/GTEMP1
                       FLHAB  = ABS(GTEMP2)
                       DELFL  = COFRM4(IS)*DT4A
                       USFM   = USNEW(IP)*MAX(FMEANWS(IP),FMEAN(IP))
                       TEMP   = USFM*DELFL
                       FLHAB  = MIN(FLHAB,TEMP)
                       IMATRAA(IS,ID,IP) = SIGN(FLHAB,GTEMP2)/DT4A*JAC
                       LIMFAC            = MIN(ONE,ABS(SIGN(FLHAB,GTEMP2/DT4A))/MAX(THR,ABS(IMATRAA(IS,ID,IP))))
                       IMATDAA(IS,ID,IP) = ZERO ! -FL(IP,ID,IS)
                     ENDIF
                   ENDDO ! ID
                 ENDDO ! IS
                 ACLOC = AC2(:,:,IP)
                 IF (.NOT. LINID) THEN
                   CALL SET_WIND( IP, WIND10, WINDTH )
                   CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
                   CALL SIN_LIN_CAV(IP,WIND10,WINDTH,FPM,IMATRA,SSINL)
                   IMATRAA(:,:,IP) = IMATRAA(:,:,IP) + SSINL
                 ENDIF
                 IF (ISHALLOW(IP) .EQ. 1) THEN
                   CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)
                   IF (MESTR .GT. 0) THEN
                     CALL TRIAD_ELDEBERKY(IP, HS, SME01, ACLOC, IMATRA, IMATDA, SSNL3, DSSNL3)
                     DO IS = 1, MSC
                       DO ID = 1, MDC
                         NEWDAC = SSNL3(IS,ID)*DT4A/MAX((1.-DT4A*DSSNL3(IS,ID)),1.)
                         MAXDAC = 0.0081*LIMFAK/(TWO*SPSIG(IS)*WK(IS,IP)**3*CG(IS,IP))*100
                         LIMFAC = ONE/MAX(ONE,NEWDAC/MAXDAC)
                         SC = SIGN(MIN(ABS(NEWDAC),MAXDAC),NEWDAC)/DT4A
                         !SSNL3(IS,ID)  = SC
                         !DSSNL3(IS,ID) = DSSNL3(IS,ID)*LIMFAC
                         !IF (ABS(SC) .GT. THR) WRITE(*,'(2I10,5F20.8)') IS, ID, NEWDAC, MAXDAC, DSSNL3(IS,ID), LIMFAC
                       END DO
                     END DO
                   ENDIF ! MESTR
                   IF (MESBR .GT. 0) CALL SDS_SWB(IP, SME01, KMWAM, ETOT, HS, ACLOC, IMATRA, IMATDA, SSBR, DSSBR)
                   IF (MESBF .GT. 0) CALL SDS_BOTF(IP,ACLOC,IMATRA,IMATDA,SSBF,DSSBF)
                   IMATDAA(:,:,IP) = IMATDAA(:,:,IP) + DSSBR  + DSSNL3 + DSSBF + DSSVEG
                   IMATRAA(:,:,IP) = IMATRAA(:,:,IP) + SSBR + SSNL3 + SSVEG
                 ENDIF ! ISHALLOW(IP) .EQ. 1
               ENDIF
             ENDIF
           ENDIF
           !IF (IP==TESTNODE) WRITE(*,*) IP, SUM(IMATRAA(:,:,IP)), DEP(IP), IOBP(IP), SUM(SSBR), SUM(SSBF), SUM(SSINL), SUM(SSNL3)
         ENDDO
!$OMP END PARALLEL
         !WRITE(*,'(A20,3F15.10)') 'FROM SPECINT', SUM(IMATRAA), SUM(IMATDAA), SUM(AC2)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SOURCE_INT_IMP_WAM_POST
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER        :: IP, IS, ID
         REAL(rkind)    :: ACLOC(MSC,MDC), VEC2RAD
         REAL(rkind)    :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         IMATDAA = 0.
         IMATRAA = 0.

         DO IP = 1, MNP
           IF ((ABS(IOBP(IP)) .NE. 1 .AND. IOBP(IP) .NE. 3)) THEN
             IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
               THWOLD(IP,1) = THWNEW(IP)
               THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))
               U10NEW(IP) = MAX(TWO,SQRT(WINDXY(IP,1)**2+WINDXY(IP,2)**2)) * WINDFAC
               Z0NEW(IP) = Z0OLD(IP,1)
               DO IS = 1, MSC
                 DO ID = 1, MDC
                   FL3(IP,ID,IS) =  AC2(IS,ID,IP) * PI2 * SPSIG(IS)
                   FL(IP,ID,IS)  =  IMATDAA(IS,ID,IP)
                   SL(IP,ID,IS)  =  IMATRAA(IS,ID,IP) * PI2 * SPSIG(IS)
                 END DO
               END DO
               CALL POSTINTRHS (FL3(IP,:,:), FL(IP,:,:), IP, IP, 1, &
     &                          THWOLD(IP,1), USOLD(IP,1), &
     &                          TAUW(IP), Z0OLD(IP,1), &
     &                          ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                          U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                          Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                          SL(IP,:,:), FCONST(IP,:), FMEANWS(IP), MIJ(IP))
               DO IS = 1, MSC
                 DO ID = 1, MDC
                   AC2(IS,ID,IP) = FL3(IP,ID,IS) / PI2 / SPSIG(IS)
                 END DO
               END DO
             END IF !
           ELSE
             IF (LSOUBOUND) THEN ! Source terms on boundary ...
               IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
                 FL = FL3
                 THWOLD(:,1) = THWNEW
                 U10NEW = MAX(TWO,SQRT(WINDXY(:,1)**2+WINDXY(:,2)**2)) * WINDFAC
                 Z0NEW(IP) = Z0OLD(IP,1)
                 THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))
                 DO IS = 1, MSC
                   DO ID = 1, MDC
                     FL3(IP,ID,IS) =  AC2(IS,ID,IP) * PI2 * SPSIG(IS)
                     FL(IP,ID,IS)  =  IMATDAA(IS,ID,IP)
                     SL(IP,ID,IS)  =  IMATRAA(IS,ID,IP) * PI2 * SPSIG(IS)
                   END DO
                 END DO
                 CALL POSTINTRHS (FL3(IP,:,:), FL(IP,:,:), IP, IP, 1, &
     &                            THWOLD(IP,1), USOLD(IP,1), &
     &                            TAUW(IP), Z0OLD(IP,1), &
     &                            ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                            U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                            Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                            SL(IP,:,:), FCONST(IP,:), FMEANWS(IP), MIJ(IP))
                 DO IS = 1, MSC
                   DO ID = 1, MDC
                     AC2(IS,ID,IP) = FL3(IP,ID,IS) / PI2 / SPSIG(IS)
                   END DO
                 END DO
               ENDIF
             ENDIF
           ENDIF
         END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INT_IP_STAT(IP,DT,LIMITER,ACLOC,ISELECT)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP, ISELECT
         LOGICAL, INTENT(IN) :: LIMITER
         REAL(rkind), INTENT(IN) :: DT
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)
         REAL(rkind)                :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER :: IS, ID

         REAL(rkind)    :: NPF(MSC)
         REAL(rkind)    :: OLDAC
         REAL(rkind)    :: NEWDAC
         REAL(rkind)    :: MAXDAC, CONST, SND, USTAR

         CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE., ISELECT, 'INT_IP_STAT')  ! 1. CALL

         CONST = PI2**2*3.0*1.0E-7*DT*SPSIG(MSC_HF(IP))
         SND   = PI2*5.6*1.0E-3

         IF (LIMITER) THEN
           IF (MELIM == 1) THEN
             NPF = 0.0081_rkind*LIMFAK/(TWO*SPSIG*WK(:,IP)**3*CG(:,IP))
           ELSE IF (MELIM == 2) THEN
             DO IS = 1, MSC
               USTAR = MAX(UFRIC(IP), G9*SND/SPSIG(IS))
               NPF(IS) = ABS((CONST*USTAR)/(SPSIG(IS)**3*WK(IS,IP)))
             END DO
           END IF
         END IF

!         if (SUM(ACLOC) .NE. SUM(ACLOC)) STOP 'NAN l. 174 wwm_specint.F90'

         DO IS = 1, MSC
           MAXDAC = NPF(IS)
           DO ID = 1, MDC
             OLDAC  = ACLOC(IS,ID)
             NEWDAC = IMATRA(IS,ID) * DT / ( 1.0 - DT * MIN(ZERO,IMATDA(IS,ID))) 
             IF (LIMITER) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
             ACLOC(IS,ID) = MAX( 0.0_rkind, OLDAC + NEWDAC ) 
           END DO
         END DO

!         write(*,'(A10,2I10,L10,I10,2F15.6)') 'after', ip, iobp(ip), limiter, iselect, sum(acloc), sum(imatra)

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE RKS_SP3(IP,DTSII,LIMITER,ACLOC,ISELECT)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP, ISELECT
         LOGICAL, INTENT(IN) :: LIMITER
         REAL(rkind), INTENT(IN)    :: DTSII
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)

         INTEGER :: IS, ID
         REAL(rkind)    :: ACOLD(MSC,MDC)
         REAL(rkind)    :: NPF(MSC)
         REAL(rkind)    :: NEWDAC, MAXDAC, CONST, SND, USTAR
         REAL(rkind)    :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         CONST = PI2**2*3.0*1.0E-7*DTSII*SPSIG(MSC_HF(IP))
         SND   = PI2*5.6*1.0E-3

         IF (LIMITER) THEN
           IF (MELIM == 1) THEN
             NPF = 0.0081_rkind*LIMFAK/(TWO*SPSIG*WK(:,IP)**3*CG(:,IP))
           ELSE IF (MELIM == 2) THEN
             DO IS = 1, MSC
               USTAR = MAX(UFRIC(IP), G9*SND/SPSIG(IS))
               NPF(IS) = ABS((CONST*USTAR)/(SPSIG(IS)**3*WK(IS,IP)))
             END DO
           END IF
         END IF
       
         CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE., ISELECT, 'RKS_SP3 1')  ! 1. CALL

         ACOLD = ACLOC

         DO IS = 1, MSC
           MAXDAC = NPF(IS)
           DO ID = 1, MDC
             NEWDAC = IMATRA(IS,ID) * DTSII / ( ONE - DTSII * MIN(ZERO,IMATDA(IS,ID)) )
             IF (LIMITER) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
             ACLOC(IS,ID) = MAX( ZERO, ACLOC(IS,ID) + NEWDAC )
           END DO
         END DO

         !WRITE(*,*) '1 RK-TVD', SUM(ACOLD), SUM(ACLOC)

         CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE., ISELECT, 'RKS_SP3 2') ! 2. CALL

         DO IS = 1, MSC
           MAXDAC = NPF(IS)
           DO ID = 1, MDC
             NEWDAC = IMATRA(IS,ID) * DTSII / ( ONE - DTSII * MIN(ZERO,IMATDA(IS,ID)) )
             IF (LIMITER) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC) 
             ACLOC(IS,ID) = MAX( ZERO, 0.75_rkind * ACOLD(IS,ID) +  0.25_rkind * ACLOC(IS,ID) + 0.25_rkind * NEWDAC)
           END DO
         END DO

         !WRITE(*,*) '2 RK-TVD', SUM(ACOLD), SUM(ACLOC)

         CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE., ISELECT, 'RKS_SP3 3') ! 3. CALL

         DO IS = 1, MSC
           MAXDAC = NPF(IS)
           DO ID = 1, MDC
             NEWDAC = IMATRA(IS,ID) * DTSII / ( ONE - DTSII * MIN(ZERO,IMATDA(IS,ID)) )
             IF (LIMITER) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
             ACLOC(IS,ID) = MAX( ZERO, ONE/THREE * ACOLD(IS,ID) + TWO/THREE * ACLOC(IS,ID) + TWO/THREE * NEWDAC)
           END DO
         END DO

         !IF (IP == 1701) WRITE(*,*) SUM(IMATRA), SUM(IMATDA), SUM(ACLOC)

         !WRITE(*,*) '3 RK-TVD', SUM(ACOLD), SUM(ACLOC)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INT_IP_DYN(IP, DT, LIMIT, DTMIN, ITRMX, ACLOC, ITER, ISELECT)

         USE DATAPOOL
#ifdef ST41
         USE W3SRC4MD_OLD
#elif ST42
         USE W3SRC4MD
#endif
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IP, ITRMX, ISELECT
         INTEGER, INTENT(OUT)       :: ITER
         LOGICAL,INTENT(IN)         :: LIMIT
         REAL(rkind), INTENT(IN)    :: DTMIN, DT
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)
         REAL(rkind)                :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER :: IS, ID

         REAL(rkind)    :: ACOLD(MSC,MDC)
         REAL(rkind)    :: NPF(MSC)
         REAL(rkind)    :: AFILT
         REAL(rkind)    :: NEWDAC, CONST, SND, DAMAX, AFAC
         REAL(rkind)    :: MAXDAC, DTMAX, DTTOT, DTLEFT, DT4SI, USTAR

         REAL(rkind), PARAMETER :: MAXDTFAC = VERYLARGE 

         LOGICAL :: LSTABLE

         REAL(rkind) :: XPP, XRR, XFILT, XREL, XFLT, FACP, DAM(MSC)

#ifdef ST_DEF
         XPP    = 0.15
         XRR    = 0.10
         XFILT  = 0.0001 ! AR: check why it must be so small ..
         XPP    = MAX ( 1.E-6_rkind , XPP )
         XRR    = MAX ( 1.E-6_rkind , XRR )
         XREL   = XRR
         XFLT   = MAX ( ZERO , XFILT ) 
         FACP   = 2*XPP / PI2 * 0.62E-3 * PI2**4 / G9**2  ! s4/m4
         DAM    = FACP / ( SPSIG * WK(:,IP)**3 ) / CG(:,IP) ! s * m³ * s4/m4 = 
         AFILT  = MAX ( DAM(MSC) , XFLT*MAXVAL(ACLOC))!*PI2*SPSIG(MSC_HF(IP)) )
#endif

         CONST = PI2**2*3.0*1.0E-7*DT*SPSIG(MSC_HF(IP))
         SND   = PI2*5.6*1.0E-3
         IMATRA = ZERO 
         IMATDA = ZERO

         IF (LIMIT) THEN
           IF (MELIM == 1) THEN
             NPF = 0.0081_rkind*LIMFAK/(TWO*SPSIG*WK(:,IP)**3*CG(:,IP))
           ELSE IF (MELIM == 2) THEN
             DO IS = 1, MSC
               USTAR = MAX(UFRIC(IP), G9*SND/SPSIG(IS))
               NPF(IS) = ABS((CONST*USTAR)/(SPSIG(IS)**3*WK(IS,IP)))
             END DO
           END IF
         END IF

         DTTOT = 0.
         ITER  = 0

         DO WHILE ( DTTOT < DT )

           ACOLD = ACLOC
           IF (ITER == 0) DT4SI = DT
           IF (ITER == 0) DTLEFT = DT
           ITER  = ITER + 1
           DTMAX =  DT

           CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE.,ISELECT, 'INT_IP_DYN')  ! 1. CALL
        
           DO ID = 1, MDC
             DO IS = 1, MSC_HF(IP)
               IF (ABS(IMATRA(IS,ID)) .GT. VERYSMALL) THEN                
                 DAMAX  = MIN ( DAM(IS) , MAX ( XREL*ACLOC(IS,ID), AFILT ) )
                 AFAC  = ONE / MAX( 1.E-10_rkind , ABS(IMATRA(IS,ID)/DAMAX) )
                 DTMAX = MIN (DTMAX ,AFAC/(MAX(1.E-10_rkind, ONE + AFAC*MIN(ZERO,IMATDA(IS,ID))))) 
               END IF
             END DO
           END DO

           DT4SI  = DTMAX
           DT4SI  = MAX(DTMIN,DT4SI) ! This is not entirely stable !!!
           DTLEFT = DT - DTTOT

           IF ( DTLEFT > THR .AND. DTLEFT < DT4SI) THEN
             DT4SI = (DT - DTTOT)
           ELSE IF ( DTLEFT .GE. DT4SI .AND. ITER .EQ. ITRMX) THEN
             DT4SI = DTLEFT 
             LSTABLE = .FALSE. 
           END IF 

           DTTOT = DTTOT + DT4SI

           IF ( DT4SI .LT. DTMAX + SMALL) THEN
             LSTABLE = .TRUE. 
           ELSE
             LSTABLE = .FALSE.
           END IF

           IF (LSTABLE) THEN
             DO IS = 1, MSC_HF(IP)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(ZERO,IMATDA(IS,ID)) )
                 ACLOC(IS,ID) = MAX( ZERO, ACLOC(IS,ID) + NEWDAC )
               END DO
             END DO
             CALL SOURCETERMS(IP,ACLOC,IMATRA,IMATDA,.FALSE.,ISELECT, 'INT_IP_DYN - STABLE1') ! 2. CALL
             DO IS = 1, MSC_HF(IP)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(ZERO,IMATDA(IS,ID)) )
                 ACLOC(IS,ID) = MAX( ZERO, 3._rkind/4._rkind * ACOLD(IS,ID) + 1._rkind/4._rkind * ACLOC(IS,ID) + 1._rkind/4._rkind * NEWDAC)
               END DO
             END DO
             CALL SOURCETERMS(IP,ACLOC,IMATRA,IMATDA, .FALSE.,ISELECT, 'INT_IP_DYN STABLE2') ! 3. CALL
             DO IS = 1, MSC_HF(IP)
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(ZERO,IMATDA(IS,ID)) )
                 ACLOC(IS,ID) = MAX( ZERO,  1._rkind/3._rkind * ACOLD(IS,ID) + 2._rkind/3._rkind * ACLOC(IS,ID) + 2._rkind/3._rkind * NEWDAC)
               END DO
             END DO
           ELSE ! .NOT. LSTABLE
             DO IS = 1, MSC_HF(IP)
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(ZERO,IMATDA(IS,ID)) )
                 IF (LIMIT) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
                 ACLOC(IS,ID) = MAX( ZERO, ACLOC(IS,ID) + NEWDAC )
               END DO
             END DO
             CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE., ISELECT, 'INT_IP_DYN - LIM1 ') ! 2. CALL
             DO IS = 1, MSC_HF(IP)
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(ZERO,IMATDA(IS,ID)) )
                 IF (LIMIT) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
                 ACLOC(IS,ID) = MAX( ZERO, 3./4. * ACOLD(IS,ID) +  1./4. * ACLOC(IS,ID) + 1./4. * NEWDAC)
               END DO
             END DO
             CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE., ISELECT, 'INT_IP_DYN - LIM2 ') ! 3. CALL
             DO IS = 1, MSC_HF(IP)
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(ZERO,IMATDA(IS,ID)) )
                 IF (LIMIT) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
                 ACLOC(IS,ID) = MAX( ZERO, 1./3. * ACOLD(IS,ID) +  2./3. * ACLOC(IS,ID) + 2./3. * NEWDAC)
               END DO
             END DO
           END IF
         END DO
         END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
         SUBROUTINE ACTION_LIMITER
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER                 :: IP, IS, ID
         REAL(rkind)             :: NEWDAC, OLDAC, NEWAC, DELT, XIMP, DELFL(MSC)
         REAL(rkind)             :: MAXDAC, CONST, SND, UFR_LIM, DELT5, USFM
         REAL(rkind)             :: ACLOC(MSC,MDC), ACOLD(MSC,MDC), MAXDACOLD


         CONST = PI2**2*3.0*1.0E-7*DT4S*SPSIG(MSC)
         SND   = PI2*5.6*1.0E-3

         DELT = DT4S
         XIMP = 1._rkind
         DELT5 = XIMP*DELT
         DELFL= COFRM4*DELT
         MAXDAC = ZERO

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(IP,ACOLD,ACLOC, MAXDAC, UFR_LIM, NEWAC, OLDAC, NEWDAC)
         DO IP = 1, MNP
           IF (DEP(IP) .LT. DMIN .OR. IOBP(IP) .EQ. 2) CYCLE
           ACOLD = AC1(:,:,IP)
           ACLOC = AC2(:,:,IP)
           DO IS = 1, MSC
             IF (MELIM .EQ. 1) THEN
               MAXDAC = 0.0081*LIMFAK/(TWO*SPSIG(IS)*WK(IS,IP)**3*CG(IS,IP))
             ELSE IF (MELIM .EQ. 2) THEN
               UFR_LIM = MAX(UFRIC(IP),G9*SND/SPSIG(IS))
               MAXDAC  = LIMFAK*ABS((CONST*UFR_LIM)/(SPSIG(IS)**3*WK(IS,IP)))
             ELSE IF (MELIM .EQ. 3) THEN
               IF (USNEW(IP) .GT. SMALL) THEN
                 MAXDAC = COFRM4(IS)*USNEW(IP)*MAX(FMEANWS(IP),FMEAN(IP))/PI2/SPSIG(IS)*DT4A
               ELSE
                 MAXDAC = 0.0081*LIMFAK/(TWO*SPSIG(IS)*WK(IS,IP)**3*CG(IS,IP))
               ENDIF
             END IF
             DO ID = 1, MDC
               NEWAC  = ACLOC(IS,ID)
               OLDAC  = ACOLD(IS,ID)
               NEWDAC = NEWAC - OLDAC
!               IF (NEWDAC .GT. 0.) THEN
                 NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)), NEWDAC)
!               ELSE
!                 IF (QBLOCAL(IP) .LT. THR) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)), NEWDAC)
!               END IF
               AC2(IS,ID,IP) = MAX( zero, OLDAC + NEWDAC )
               !IF (MELIM .EQ. 3) FL3(IP,ID,IS) = AC2(IP,IS,ID) * PI2 * SPSIG(IS)
             END DO
           END DO
         END DO
         END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
         SUBROUTINE ACTION_LIMITER_LOCAL(IP,ACLOC,ACOLD)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IP
         INTEGER                    :: IS, ID
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(IN)    :: ACOLD(MSC,MDC)
         REAL(rkind)                :: NEWDAC, OLDAC, NEWAC, DELT, XIMP, DELFL(MSC)
         REAL(rkind)                :: MAXDAC, CONST, SND, UFR_LIM, DELT5, USFM
         REAL(rkind)                :: MAXDACOLD

         CONST = PI2**2*3.0*1.0E-7*DT4S*SPSIG(MSC)
         SND   = PI2*5.6*1.0E-3

         DELT = DT4S
         XIMP = 1._rkind
         DELT5 = XIMP*DELT
         DELFL= COFRM4*DELT
         MAXDAC = ZERO

         DO IS = 1, MSC

           IF (MELIM .EQ. 1) THEN
             MAXDAC = 0.0081*LIMFAK/(TWO*SPSIG(IS)*WK(IS,IP)**3*CG(IS,IP))
           ELSE IF (MELIM .EQ. 2) THEN
             UFR_LIM = MAX(UFRIC(IP),G9*SND/SPSIG(IS))
             MAXDAC  = LIMFAK*ABS((CONST*UFR_LIM)/(SPSIG(IS)**3*WK(IS,IP)))
           ELSE IF (MELIM .EQ. 3) THEN
             IF (USNEW(IP) .GT. SMALL) THEN
               MAXDAC = COFRM4(IS)*USNEW(IP)*MAX(FMEANWS(IP),FMEAN(IP))/PI2/SPSIG(IS)*DT4A
             ELSE
               MAXDAC = 0.0081*LIMFAK/(TWO*SPSIG(IS)*WK(IS,IP)**3*CG(IS,IP))
             ENDIF
           END IF

           DO ID = 1, MDC
             NEWAC  = ACLOC(IS,ID)
             OLDAC  = ACOLD(IS,ID)
             NEWDAC = NEWAC - OLDAC
             NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)), NEWDAC)
             ACLOC(IS,ID) = MAX( zero, OLDAC + NEWDAC )
           END DO

         END DO

         END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
         SUBROUTINE BREAK_LIMIT(IP,ACLOC,SSBRL)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)  :: IP

         REAL(rkind), INTENT(INOUT)  :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(OUT)    :: SSBRL(MSC,MDC)

         REAL(rkind)                 :: EFTAIL, HS
         REAL(rkind)                 :: EMAX, RATIO, ETOT
         REAL(rkind)                 :: DINTSPEC

         ETOT   = 0.0
         EFTAIL = 1.0 / (PTAIL(1)-1.0)

         ETOT = DINTSPEC(IP,ACLOC)

         HS = 4.*SQRT(ETOT)

         !IF (LMONO_IN) HMAX(IP) = HMAX(IP) * SQRT(2.)

         EMAX = 1./16. * (HMAX(IP))**2

         IF (ETOT .GT. EMAX .AND. ETOT .GT. THR) THEN
           RATIO = EMAX/ETOT
           SSBRL = ACLOC - RATIO * ACLOC
           ACLOC = RATIO * ACLOC
         END IF
         END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
         SUBROUTINE BREAK_LIMIT_ALL
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER              :: IP
         REAL(rkind)          :: HS
         REAL(rkind)          :: EMAX, RATIO, ETOT
         REAL(rkind)          :: DINTSPEC
         REAL(rkind)          :: ACLOC(MSC, MDC)
!      Print *, 'Passing BREAK_LIMIT_ALL'
         DO IP = 1, MNP
           ACLOC = AC2(:,:,IP)
           IF (ISHALLOW(IP) .EQ. 0) CYCLE
           ETOT = DINTSPEC(IP,ACLOC)
           HS = 4.*SQRT(ETOT)
           EMAX = 1./16. * (HMAX(IP))**2
!        WRITE(300,*) 'IP=', IP, ' HMAX=', HMAX(IP), ' DEP=', DEP(IP)
!        WRITE(300,*) '   ', IP, ' EMAX=', EMAX, ' ETOT=', ETOT
!        WRITE(300,*) '   ', IP, ' HS=', HS, ' BRHD=', BRHD

           IF (ETOT .GT. EMAX) THEN
             if(myrank==0) WRITE(300,*) '   break XP=', XP(IP)
             RATIO = EMAX/ETOT
             AC2(:,:,IP) = RATIO * ACLOC(:,:)
             AC1(:,:,IP) = RATIO * ACLOC(:,:)
           END IF
         END DO
         END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
         SUBROUTINE RESCALE_SPECTRUM()
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP, IS, ID
         REAL(rkind)    :: ETOTAL, EPOSIT
         REAL(rkind)    :: FACTOR
         LOGICAL :: ENEG

         DO IP = 1, MNP
            DO IS = 1, MSC
               ETOTAL = ZERO 
               EPOSIT = ZERO 
               ENEG   = .FALSE.
               DO ID = 1, MDC
                  ETOTAL = ETOTAL + AC2(IS,ID,IP)
                  IF (AC2(IS,ID,IP) > ZERO) THEN
                     EPOSIT = EPOSIT + AC2(IS,ID,IP)
                  ELSE
                     ENEG = .TRUE.
                  END IF
               END DO
               IF (ENEG) THEN
                  IF (EPOSIT .GT. VERYSMALL) THEN
                    FACTOR = ETOTAL/EPOSIT
                  ELSE 
                    FACTOR = 0.
                  END IF
                  DO ID = 1, MDC
                     IF (AC2(IS,ID,IP) < ZERO) AC2(IS,ID,IP) = ZERO 
                     IF (FACTOR >= ZERO)  AC2(IS,ID,IP) = AC2(IS,ID,IP)*FACTOR
                     AC2(IS,ID,IP) = MAX(zero,AC2(IS,ID,IP))
                  END DO
               END IF
            END DO
         END DO
         END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
         SUBROUTINE SETSHALLOW
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP

         DO IP = 1, MNP
           IF (WK(1,IP)*DEP(IP) .LT. PI) THEN
             ISHALLOW(IP) = 1
           ELSE
             ISHALLOW(IP) = 0
           END IF
         END DO
         END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
