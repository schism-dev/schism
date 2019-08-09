      PROGRAM STATISTICS
      USE STAT_POOL 
      IMPLICIT NONE
      
      REAL*8, ALLOCATABLE    :: ACLOC(:,:), ACLOC2D(:,:)
      REAL*8, ALLOCATABLE    :: WKLOC(:)
      REAL*8, ALLOCATABLE  :: TMPE(:), TMPES(:)

      REAL*8               :: DEPLOC, CURTXYLOC(2), WINDXY(2)
      REAL*8               :: HS,TM01,TM02,KLM,WLM,ETMPO,ETMPS

      REAL*8                :: ETOT, DELTAE, ETOT1

      CHARACTER(LEN = 30)  :: CTMP30 
      CHARACTER(LEN = 15)  :: CTMP15
      CHARACTER(LEN = 5)   :: NAMEBUOY 

      INTEGER              :: IS, ID, IB, I, J, K, IT, ISMAX, IBUOYS
      INTEGER              :: IS_O, IT_O
      INTEGER              :: dateFMT = 0

      INTEGER, ALLOCATABLE :: I2DSPECOUT(:)

      INTEGER, ALLOCATABLE :: IFIND(:), ICOUNT(:), FILETYPE(:)
      
      PI  = 4.* DATAN(1.d0)
      PI2 = 2.* PI

!------------------------------------------------------------------------      
! open files  
!------------------------------------------------------------------------
      
      CALL TEST_FILE_EXIST_DIE("Missing station list : ", "ndbcdaten.dat")
      OPEN(1,  FILE = 'ndbcdaten.dat', STATUS='OLD') 
      CALL TEST_FILE_EXIST_DIE("Missing list of spectra : ", "stationen.dat")
      OPEN(2,  FILE = 'stationen.dat', STATUS='OLD')
      CALL TEST_FILE_EXIST_DIE("Missing list of setup : ", "setup.dat")
      OPEN(3,  FILE = 'setup.dat',     STATUS='OLD')

!AR: todo get rid of it and read directly from the buoys ...

!      WRITE(*,*) 'Please define desirted cut-off freq. for integration'
!      READ(*,*) CUT_OFF 
!      WRITE(*,*) 'Please define ramptime for output'
!      READ(*,*) RAMPTIME 
      
      WRITE(*,*) 'Please define date format'
      WRITE(*,*) '(1) YYYY MM DD hh'
      WRITE(*,*) '(2) YYYY MM DD hh mm'
      CUT_OFF = 0.4 
      dateFMT = 1. 
!      READ(*,*) dateFMT
      
      BUOYS = 0  

!------------------------------------------------------------------------      
! get sizes 
!------------------------------------------------------------------------
      
      DO
        READ (2, *,IOSTAT = DATAEND) 
        IF (DATAEND < 0) THEN
          EXIT
        END IF
        BUOYS = BUOYS  + 1
      END DO 
      WRITE(1001,*) BUOYS
      REWIND 2
      ALLOCATE (STATIONNAMES(BUOYS), IFIND(BUOYS), ICOUNT(BUOYS), FILETYPE(BUOYS), I2DSPECOUT(BUOYS))

      I2DSPECOUT(:) = 0
         
      DO IB = 1, BUOYS  
        READ (2,'(A15)',IOSTAT = DATAEND) STATIONNAMES(IB) 
        WRITE(1001,'(A15)',IOSTAT = DATAEND) STATIONNAMES(IB)
      END DO  
 
      ALLOCATE (BNAMES(BUOYS)) 
      
      DO IB = 1, BUOYS 
        CTMP30 = TRIM(STATIONNAMES(IB))
        BNAMES(IB) = CTMP30(1:5) 
        WRITE(1002,*)  IB, BNAMES(IB), TRIM(STATIONNAMES(IB))
      END DO

      CALL TEST_FILE_EXIST_DIE("Missing station names : ", TRIM(STATIONNAMES(1)))
      OPEN(4, FILE = TRIM(STATIONNAMES(1)), STATUS='OLD', FORM='UNFORMATTED')

      READ (4) MSC
      READ (4) MDC
      ALLOCATE (SPSIG(MSC),SPDIR(MDC),ACLOC(MSC,MDC),ACLOC2D(MSC,MDC),DFREQ_SP(MSC),COSTH(MDC),SINTH(MDC),FREQ_SP(MSC))
      ALLOCATE (AUX_M0_SP(MSC), AUX_M1_SP(MSC), AUX_M2_SP(MSC), E_SP(MSC))
      ALLOCATE (DS_INCR_SP(0:MSC+1))
      ALLOCATE (TMPES(MSC)); TMPES = 0.d0
      READ (4) SPSIG
      READ (4) SPDIR
      READ (4) IFIND(1)
      DDIR = SPDIR(2)-SPDIR(1)

      DO ID = 1, MDC
        COSTH(ID) = COS(SPDIR(ID))
        SINTH(ID) = SIN(SPDIR(ID))
      END DO

      FREQ_SP = SPSIG / PI2

      WRITE(2000,*) FREQ_SP

      DO IS = 1, MSC - 1
        DFREQ_SP(IS) = FREQ_SP(IS+1) - FREQ_SP(IS) 
      END DO

      WRITE(2000,*) DFREQ_SP

      SGLOW  = SPSIG(1)
      SGHIG  = SPSIG(MSC)
      FRINTF = DLOG(SGHIG/SGLOW)/REAL(MSC-1)

      ISMAX = MSC 
      DO IS = 1, MSC
        IF (FREQ_SP(IS) .GT. CUT_OFF) THEN
          ISMAX = IS 
        EXIT
        END IF
      END DO

      WRITE(2000,*) 'ISMAX', ISMAX

!------------------------------------------------------------------------      
! number of time steps in simulations ... 
!------------------------------------------------------------------------

      N_DT_SP = 0
      DO
        READ (4,IOSTAT = DATAEND) CTMP15
        READ (4,IOSTAT = DATAEND) DEPLOC
        READ (4,IOSTAT = DATAEND) CURTXYLOC
        READ (4,IOSTAT = DATAEND) ACLOC
        READ (4,IOSTAT = DATAEND) ACLOC2D
        IF (DATAEND < 0) THEN
          EXIT
        END IF
        N_DT_SP = N_DT_SP + 1
      END DO 

      REWIND (4)
      CLOSE(4)  

!------------------------------------------------------------------------      
! alloc and init sim arrays ... 
!------------------------------------------------------------------------

      ALLOCATE (DATUM_SP (N_DT_SP))
      ALLOCATE (ZEIT_SP  (N_DT_SP))
      ALLOCATE (HS_SP    (N_DT_SP, BUOYS))
      ALLOCATE (TM01_SP  (N_DT_SP, BUOYS))
      ALLOCATE (TM02_SP  (N_DT_SP, BUOYS))
      ALLOCATE (EMAX_SP  (N_DT_SP, BUOYS))
      ALLOCATE (TP_SP    (N_DT_SP, BUOYS))

      DATUM_SP = ''
      ZEIT_SP  = 0.d0
      HS_SP    = 0.
      TM01_SP  = 0. 
      TM02_SP  = 0.

!------------------------------------------------------------------------      
! read simulation results ... 
!------------------------------------------------------------------------

      IBUOYS = 0
      DO IB = 1, BUOYS
        CALL TEST_FILE_EXIST_DIE("Missing station file : ", TRIM(STATIONNAMES(IB)))
        OPEN(4, FILE = TRIM(STATIONNAMES(IB)), STATUS='OLD', FORM='UNFORMATTED')
        READ(4) MSC
        READ(4) MDC
        READ(4) SPSIG
        READ(4) SPDIR
        READ(4) IFIND(IB)
        WRITE(2220,*) TRIM(STATIONNAMES(IB))
        WRITE(2221,*) TRIM(STATIONNAMES(IB))
        WRITE(1003,*) 'HEADER INFORMATIONS OF THE CERTAIN BUOY'
        WRITE(1003,*) IB, IFIND(IB), TRIM(STATIONNAMES(IB))
        WRITE(1003,*) MSC, MDC
        WRITE(1003,*) SPSIG, SPDIR
        WRITE(1003,*) 'TIME HISTORY'

        IF (IB == 1) THEN
          WRITE(4002,*) '2D SPECTRA WWMIII'
          WRITE(4002,*) 'START WRITING 2D SPECTRA FOR STATION =', IB, TRIM(STATIONNAMES(IB))
          WRITE(4002,*) 'NUMBER OF FREQ', MSC
          WRITE(4002,*) 'NUMBER OF DIRECTIONS', MDC
          WRITE(4002,*) 'SPECTRAL GRID IN RAD AND HZ' 
          DO IS = 1, MSC
            WRITE(4002,*) IS, SPSIG(IS), FREQ_SP(IS)
          END DO 
          WRITE(4002,*) 'DIRECTIONAL GRID DDIR =',DDIR
          DO ID = 1, MDC
            WRITE(4002,*) ID, SPDIR(ID)
          END DO
        ENDIF
        IF (IFIND(IB) .NE. 0) THEN
          IBUOYS = IBUOYS + 1
          DO IT = 1, N_DT_SP
            READ(4) DATUM_SP(IT)
            READ(4) DEPLOC
            READ(4) CURTXYLOC
            READ(4) ACLOC
            READ(4) ACLOC2D
            WRITE(2000,*) IB, IT, 'CURTXYLOC', CURTXYLOC
            CALL CT2MJD(DATUM_SP(IT), XMJD)
            WRITE(2001,*) DATUM_SP(IT), XMJD
            ZEIT_SP(IT)  = XMJD
            AUX_M0_SP(:) = 0.0
            AUX_M1_SP(:) = 0.0
            AUX_M2_SP(:) = 0.0
            E_SP = 0.0
            DO K = 1, ISMAX
              E_SP(K) = SUM(ACLOC2D(K,:)) * DDIR * PI2 !* SPSIG(K) 
            END DO
            ETOT1  = SUM(SUM(ACLOC(:,:),DIM=2) * SPSIG(:)**2 * FRINTF*DDIR )
            IF (IB == 1) THEN
              WRITE(4002,*) IB, IT, TRIM(DATUM_SP(IT)), ZEIT_SP(IT)
              DO IS = 1, MSC
                DO ID = 1, MDC
                  write(4002,*) FREQ_SP(IS), SPDIR(ID), ACLOC2D(IS,ID)*PI2
                ENDDO
              ENDDO
            ENDIF
            EMAX_SP(IT,IB)  = E_SP(1)
            ISEMAX          = 1
            DO K = 1, ISMAX - 1
              DELTAE = 0.5*(E_SP(K+1)+E_SP(K))
              AUX_M0_SP(K) = DELTAE * DFREQ_SP(K)
              AUX_M1_SP(K) = DELTAE * DFREQ_SP(K) * FREQ_SP(K)
              AUX_M2_SP(K) = DELTAE * DFREQ_SP(K) * FREQ_SP(K)**2.
              IF (E_SP(K) .GE. EMAX_SP(IT,IB)) THEN
                EMAX_SP(IT,IB) = E_SP(K)
                ISEMAX = K
              ENDIF
              !write(*,*) ISEMAX, E_SP(K), EMAX_SP(IT,IB)
            END DO
            IF (E_SP(ISMAX) .GT. EMAX_SP(IT,IB)) EMAX_SP(IT,IB) = E_SP(ISMAX)
            M0 = SUM(AUX_M0_SP(:))
            M1 = SUM(AUX_M1_SP(:))
            M2 = SUM(AUX_M2_SP(:))
            !write(*,*) etot1, m0
            WRITE(2000,*) 'M0, M1, M2', M0, M1, M2
            IF (M0 .GT. 1.E-14) THEN
              HS_SP  (IT,IB) =  4*M0**0.5
              TM01_SP(IT,IB) = (M0/M1)
              TM02_SP(IT,IB) = (M0/M2)**0.5
              TP_SP  (IT,IB) = 1./FREQ_SP(ISEMAX)
            ELSE
              HS_SP  (IT,IB) = 0.
              TM01_SP(IT,IB) = 0.
              TM02_SP(IT,IB) = 0.
              TP_SP  (IT,IB) = 0.
            END IF
          END DO
        ELSE
          HS_SP  (:,IB) = 0. 
          TM01_SP(:,IB) = 0. 
          TM02_SP(:,IB) = 0. 
          TP_SP  (:,IB) = 0.
        END IF
        CLOSE(4)
      END DO 

!------------------------------------------------------------------------      
! read observation station namelist ...
!------------------------------------------------------------------------

      BUOYS_O = 0  
      DO
        READ (1, '(A)',IOSTAT = DATAEND) DUMP 
        IF (DATAEND < 0) THEN
          EXIT
        END IF         
        BUOYS_O = BUOYS_O  + 1       
      END DO
      REWIND 1 
      
      WRITE(*,*) 'Es sind', BUOYS_O, ' Bojen vorhanden' 
      
      IF (BUOYS_O .NE. BUOYS) THEN
        WRITE(*,*) 'Take care ther are less observation points than output locations in the station list'
      END IF

!------------------------------------------------------------------------      
! alloc observation arrays ... 
!------------------------------------------------------------------------

      ALLOCATE (STATIONNAMES_O(BUOYS_O))
      ALLOCATE (FREQ_O(MSC_O_MAX, BUOYS_O))
      ALLOCATE (DFREQ_O(MSC_O_MAX, BUOYS_O)) 
      ALLOCATE (FMT_O_F(BUOYS_O))
      ALLOCATE (FMT_O(BUOYS_O))
      
      FREQ_O(:,:)  = 0.0
      DFREQ_O(:,:)  = 0.0
      STATIONNAMES_O = ''
      FMT_O_F = ''
      FMT_O = '' 

!------------------------------------------------------------------------      
! read station names ... 
!------------------------------------------------------------------------
      
      DO IB = 1, BUOYS_O
        READ (1,*) STATIONNAMES_O(IB), FILETYPE(IB)
        WRITE(1004,*) STATIONNAMES_O(IB), FILETYPE(IB)
      END DO  
      
      ALLOCATE (MSC_O(BUOYS_O))
      
      MSC_O(:) = 0
      
      DO I = 1, BUOYS_O
        READ (3,*) MSC_O(I)
        WRITE(CHTMP,999) MSC_O(I) 
        FMT_O_F(I) = '('//CHTMP//'F7.3'//')'   ! Frequency Format 
        FMT_O(I)   = '('//CHTMP//'F7.2'//')'   ! Energy Format
        WRITE(1005,*) I, MSC_O(I), BUOYS_O, FMT_O_F(I), FMT_O(I)
      END DO 
      CLOSE(3)
      
      ALLOCATE (N_DT_O(BUOYS_O))
      N_DT_O = 0

!------------------------------------------------------------------------      
! estimate number of measurments in time ... 
!------------------------------------------------------------------------

      DO IB = 1, BUOYS_O
        CALL TEST_FILE_EXIST_DIE("Missing observation : ", TRIM(STATIONNAMES_O(IB)))
        OPEN(5, FILE = TRIM(STATIONNAMES_O(IB)), STATUS='OLD')  
        IF (FILETYPE(IB) .EQ. 1) THEN
          ALLOCATE (TMPE(MSC_O(IB)))
          if(dateFMT == 1) then
            READ(5,  '(A14'//','//FMT_O_F(IB)//')') DUMP , TMPE(:) !
          else if(dateFMT == 2) then
            READ(5,  '(A17'//','//FMT_O_F(IB)//')') DUMP , TMPE(:) !
          end if
          DO IS_O = 1, MSC_O(IB)
            FREQ_O(IS_O,IB) = TMPE(IS_O)
            !write(*,*) IB, IS_O, TMPE(IS_O)
            if (TMPE(IS_O) .le. 0.) then
              write(*,*) 'pls. check setup.dat', ib, is_o
              stop
            endif
          END DO
          DEALLOCATE(TMPE)
          DO IS_O = 1, MSC_O_MAX - 1
            DFREQ_O(IS_O,IB) =  FREQ_O(IS_O+1,IB) - FREQ_O(IS_O,IB)
          END DO
          REWIND(5)
          N_DT_O(IB) = 0 
          DO
            READ (5, '(A500)',IOSTAT = DATAEND) DUMP
            IF (DATAEND < 0) THEN
              EXIT
            END IF
            N_DT_O(IB) = N_DT_O(IB) + 1
          END DO 
          N_DT_O(IB) = N_DT_O(IB) - 1
!          WRITE (*,*) 'Anzahl der Zeitschritte der Messungen', N_DT_O(IB)
          CLOSE(5)
        ELSE IF (FILETYPE(IB) .EQ. 2) THEN
          N_DT_O(IB) = 0
          DO
            READ (5, '(A500)',IOSTAT = DATAEND) DUMP
            IF (DATAEND < 0) THEN
              EXIT
            END IF
            N_DT_O(IB) = N_DT_O(IB) + 1
          END DO
          N_DT_O(IB) = N_DT_O(IB) - 1        
        ENDIF
      END DO     
      N_DT_O_MAX = MAXVAL(N_DT_O(:))

!------------------------------------------------------------------------      
! allocate arrays needed for statistics ... 
!------------------------------------------------------------------------
      
      ALLOCATE (E_O        (MSC_O_MAX)); E_O = 0.0
      ALLOCATE (E_OS       (N_DT_O_MAX, MSC_O_MAX)); E_OS = 0.0
      ALLOCATE (E_OSF      (N_DT_O_MAX, MSC)); E_OSF = 0.
      ALLOCATE (E_BT       (N_DT_O_MAX, MSC)); E_BT = 0.
      ALLOCATE (RMS_B      (BUOYS,MSC)); RMS_B = 0.
      ALLOCATE (AUX_M0_O   (MSC_O_MAX)); AUX_M0_O = 0.
      ALLOCATE (AUX_M1_O   (MSC_O_MAX)); AUX_M1_O = 0.
      ALLOCATE (AUX_M2_O   (MSC_O_MAX)); AUX_M2_O = 0.
      ALLOCATE (DATUM_O    (N_DT_O_MAX, BUOYS_O)); DATUM_O = '' 
      ALLOCATE (ZEIT_O     (N_DT_O_MAX, BUOYS_O)); ZEIT_O = 0.
      ALLOCATE (HS_O       (N_DT_O_MAX, BUOYS_O)); HS_O = 0.
      ALLOCATE (TM01_O     (N_DT_O_MAX, BUOYS_O)); TM01_O = 0.
      ALLOCATE (TM02_O     (N_DT_O_MAX, BUOYS_O)); TM02_O = 0.   
      ALLOCATE (TP_O     (N_DT_O_MAX, BUOYS_O)); TP_O = 0.
      ALLOCATE (EMAX_O     (N_DT_O_MAX, BUOYS_O)); EMAX_O = 0.

!------------------------------------------------------------------------      
! read the observed wave spectra and estimate spectral moments ...
!------------------------------------------------------------------------
            
      DO IB = 1, BUOYS_O
        IF (IFIND(IB) .EQ. 0) CYCLE
        CALL TEST_FILE_EXIST_DIE("Missing obs file : ", TRIM(STATIONNAMES_O(IB)))
        OPEN(5, FILE = TRIM(STATIONNAMES_O(IB)), STATUS='OLD')
        !WRITE(2222,*) STATIONNAMES_O(IB)!, N_DT_O(IB)
        !WRITE(2223,*) STATIONNAMES_O(IB)!, N_DT_O(IB)
         READ (5, '(A500)') DUMP 
!        WRITE(*,*) DUMP
        DO I = 1, N_DT_O(IB)

          IF (FILETYPE(IB) .EQ. 1) THEN
            ALLOCATE (TMPE(MSC_O(IB)))
            IF(dateFMT == 1) THEN
!             YYYY MM DD hh
              READ  (5, '(A4,A,A2,A,A2,A,A2'//','//FMT_O(IB)//')') YYYY,DUMP3,MM,DUMP3,DD,DUMP3,hh, TMPE(:)
              !WRITE(2222, '(A4,A,A2,A,A2,A,A2)') YYYY,DUMP3,MM,DUMP3,DD,DUMP3,hh
              DATUM_O(I,IB) = YYYY//MM//DD//'.'//HH//'0000'
            ELSE IF(dateFMT == 2) THEN
!             YYYY MM DD hh MINUTE
              READ  (5, '(A4,A,A2,A,A2,A,A2,A,A2'//','//FMT_O(IB)//')') YYYY,DUMP3,MM,DUMP3,DD,DUMP3,hh,DUMP3,MINUTE, TMPE(:)
              DATUM_O(I,IB) = YYYY//MM//DD//'.'//HH//MINUTE//'00'
            END IF
            DO IS_O = 1, MSC_O(IB)
              E_O(IS_O) = TMPE(IS_O)
            END DO
!           WRITE(*,*) MSC_O(IB), FMT_O(IB)
            DEALLOCATE(TMPE)
!           WRITE(*,*) 'ZEIT', I, N_DT_O(IB)
            !WRITE(*, '(A14'//','//FMT_O//')') DATUM_O(I,IB)!, E_O(:)          
            !write(*,*) DATUM_O(I,IB)
            CALL CT2MJD(DATUM_O(I,IB), XMJD)
            ZEIT_O(I,IB) = XMJD
            WRITE(2002,*) IB, I, DATUM_O(I,IB), ZEIT_O(I,IB) 
            AUX_M0_O(:) = 0.0
            AUX_M1_O(:) = 0.0
            AUX_M2_O(:) = 0.0
            ISEMAX = 1
            EMAX_O(I,IB) = E_O(1) 
            DO K = 1, MSC_O(IB) - 1
              IF (FREQ_O(K,IB) .LE. CUT_OFF) THEN
                IF (E_O(K) .LT. 0.0) E_O(K) = 0.0
                DELTAE = 0.5*(E_O(K+1)+E_O(K))
                AUX_M0_O(K) = DELTAE * DFREQ_O(K,IB)
                AUX_M1_O(K) = DELTAE * DFREQ_O(K,IB) * FREQ_O(K,IB)
                AUX_M2_O(K) = DELTAE * DFREQ_O(K,IB) * FREQ_O(K,IB)**2.
                IF (E_O(K) .GE. EMAX_O(I,IB)) THEN
                  EMAX_O(I,IB) = E_O(K)
                  ISEMAX = K 
                ENDIF
!               WRITE(*,*) AUX_M0_O(K), AUX_M1_O(K), AUX_M2_O(K)
              ELSE 
                AUX_M0_O(K) = SMALL 
                AUX_M1_O(K) = SMALL 
                AUX_M2_O(K) = SMALL 
              END IF
            END DO          
            IF (E_O(ISMAX) .GT. EMAX_O(I,IB)) EMAX_O(I,IB) = E_O(ISMAX)
            M0 = SUM(AUX_M0_O(:))
            M1 = SUM(AUX_M1_O(:))
            M2 = SUM(AUX_M2_O(:))
            IF (M0 .GT. 1.E-14) THEN
              HS_O  (I,IB) =  4.*(M0)**0.5
              TM01_O(I,IB) = (M0/M1)
              TM02_O(I,IB) = (M0/M2)**0.5
              TP_O(I,IB)   = 1./FREQ_O(ISEMAX,IB)
            ELSE
              HS_O  (I,IB) = 0.
              TM01_O(I,IB) = 0. 
              TM02_O(I,IB) = 0.
              TP_O(I,IB)   = 0. 
            END IF
            WRITE(1006,'(4F15.8)') HS_O  (I,IB), TM01_O(I,IB), TM02_O(I,IB), TP_O(I,IB)
          ELSE IF (FILETYPE(IB) .EQ. 2) THEN 
            IF(dateFMT == 1) THEN
!            YYYY MM DD hh
              READ  (5, '(A4,A,A2,A,A2,A,A2'//','//FMT_O(IB)//')') YYYY,DUMP3,MM,DUMP3,DD,DUMP3,hh,HS_O(I,IB),TP_O(I,IB),DUMP_R,DUMP_R,DUMP_R,DUMP_R,DUMP_R,DUMP_R,DUMP_R,DUMP_R
              !WRITE(2222, '(A4,A,A2,A,A2,A,A2)') YYYY,DUMP3,MM,DUMP3,DD,DUMP3,hh
              DATUM_O(I,IB) = YYYY//MM//DD//'.'//HH//'0000'
              CALL CT2MJD(DATUM_O(I,IB), XMJD)
              ZEIT_O(I,IB) = XMJD
              TM01_O(I,IB) = 0.
              TM02_O(I,IB) = 0.
            ELSE IF(dateFMT == 2) THEN
!             YYYY MM DD hh MINUTE
              !READ  (5, '(A4,A,A2,A,A2,A,A2,A,A2'//','//FMT_O(IB)//')') YYYY,DUMP3,MM,DUMP3,DD,DUMP3,hh,DUMP3,MINUTE,DUMP_R,DUMP_R,DUMP_R,HS_O(I,IB),TP_O(I,IB),DUMP_R,DUMP_R,DUMP_R,DUMP_R,DUMP_R,DUMP_R,DUMP_R,DUMP_R
              READ  (5, '(A4,A,A2,A,A2,A,A2,A,A2,I4,F5.1,F5.1,2F6.2,F6.2,I4,F7.1,F6.1,2F6.1,F5.1,F6.2)') YYYY,DUMP3,MM,DUMP3,DD,DUMP3,hh,DUMP3,MINUTE,DUMP_I,DUMP_R,DUMP_R,HS_O(I,IB),TP_O(I,IB),DUMP_R,DUMP_I,DUMP_R,DUMP_R,DUMP_R, DUMP_R, DUMP_R, DUMP_R
              !write(*,*) HS_O(I,IB),TP_O(I,IB)
              !2005 01 01 00 00 190 10.0 11.0  1.30  5.00 99.00 999 9999.0   7.2   5.6 999.0  1.4 99.00
              !WRITE  (*,'(A4,A,A2,A,A2,A,A2,A,A2,I4,F5.1,F5.1,2F6.2,F6.2,I4,F7.1,F6.1,2F6.1,F5.1,F6.2)') YYYY,DUMP3,MM,DUMP3,DD,DUMP3,hh,DUMP3,MINUTE,DUMP_I,DUMP_R,DUMP_R,HS_O(I,IB),TP_O(I,IB),DUMP_R,DUMP_I,DUMP_R,DUMP_R,DUMP_R, DUMP_R, DUMP_R, DUMP_R
               
              DATUM_O(I,IB) = YYYY//MM//DD//'.'//HH//MINUTE//'00'
              CALL CT2MJD(DATUM_O(I,IB), XMJD)
              ZEIT_O(I,IB) = XMJD
              TM01_O(I,IB) = 0.
              TM02_O(I,IB) = 0.
            END IF
          ENDIF
        END DO
        CLOSE(5)
!        PAUSE
      END DO

      WRITE(*,*) 'MIN MAX ZEITEN', MINVAL(ZEIT_SP), MAXVAL(ZEIT_SP)
!      
!------------------------------------------------------------------------      
! allocate more arrays for stats's 
!------------------------------------------------------------------------
         
      ALLOCATE (DIFF_HS_SP       (N_DT_O_MAX,BUOYS)); DIFF_HS_SP = 0.
      ALLOCATE (DIFF_TM01_SP     (N_DT_O_MAX,BUOYS)); DIFF_TM01_SP = 0.
      ALLOCATE (DIFF_TM02_SP     (N_DT_O_MAX,BUOYS)); DIFF_TM02_SP = 0.
      ALLOCATE (DIFF2_HS_SP      (N_DT_O_MAX,BUOYS)); DIFF2_HS_SP = 0.
      ALLOCATE (DIFF2_TM01_SP    (N_DT_O_MAX,BUOYS)); DIFF2_TM01_SP = 0.
      ALLOCATE (DIFF2_TM02_SP    (N_DT_O_MAX,BUOYS)); DIFF2_TM02_SP = 0.
      ALLOCATE (DIFF3_HS_SP      (N_DT_O_MAX,BUOYS)); DIFF3_HS_SP = 0.
      ALLOCATE (DIFF3_TM01_SP    (N_DT_O_MAX,BUOYS)); DIFF3_TM01_SP = 0.
      ALLOCATE (DIFF3_TM02_SP    (N_DT_O_MAX,BUOYS)); DIFF3_TM02_SP = 0.
      ALLOCATE (DIFF4_HS_SP      (N_DT_O_MAX,BUOYS)); DIFF4_HS_SP = 0.
      ALLOCATE (DIFF4_TM01_SP    (N_DT_O_MAX,BUOYS)); DIFF4_TM01_SP = 0.
      ALLOCATE (DIFF4_TM02_SP    (N_DT_O_MAX,BUOYS)); DIFF4_TM02_SP = 0.      
      ALLOCATE (ABS_DIFF_HS_SP   (N_DT_O_MAX,BUOYS)); ABS_DIFF_HS_SP = 0.
      ALLOCATE (ABS_DIFF_TM01_SP (N_DT_O_MAX,BUOYS)); ABS_DIFF_TM01_SP = 0.
      ALLOCATE (ABS_DIFF_TM02_SP (N_DT_O_MAX,BUOYS)); ABS_DIFF_TM02_SP = 0.
      ALLOCATE (HS_SP_CLEAN      (N_DT_O_MAX,BUOYS)); HS_SP_CLEAN = 0.
      ALLOCATE (TM01_SP_CLEAN    (N_DT_O_MAX,BUOYS)); TM01_SP_CLEAN = 0.
      ALLOCATE (TM02_SP_CLEAN    (N_DT_O_MAX,BUOYS)); TM02_SP_CLEAN = 0.
      ALLOCATE (HS_O_CLEAN       (N_DT_O_MAX,BUOYS)); HS_O_CLEAN = 0.
      ALLOCATE (TM01_O_CLEAN     (N_DT_O_MAX,BUOYS)); TM01_O_CLEAN = 0.
      ALLOCATE (TM02_O_CLEAN     (N_DT_O_MAX,BUOYS)); TM02_O_CLEAN = 0.

      ALLOCATE (DIFF_TP_SP       (N_DT_O_MAX,BUOYS)); DIFF_TP_SP = 0.
      ALLOCATE (DIFF2_TP_SP      (N_DT_O_MAX,BUOYS)); DIFF2_TP_SP = 0.
      ALLOCATE (DIFF3_TP_SP      (N_DT_O_MAX,BUOYS)); DIFF3_TP_SP = 0.
      ALLOCATE (DIFF4_TP_SP      (N_DT_O_MAX,BUOYS)); DIFF4_TP_SP = 0.
      ALLOCATE (ABS_DIFF_TP_SP   (N_DT_O_MAX,BUOYS)); ABS_DIFF_TP_SP = 0.
      ALLOCATE (TP_O_CLEAN       (N_DT_O_MAX,BUOYS)); TP_O_CLEAN = 0.
      ALLOCATE (TP_SP_CLEAN      (N_DT_O_MAX,BUOYS)); TP_SP_CLEAN = 0.


      ALLOCATE (N_OUT_TIME     (BUOYS)); N_OUT_TIME = 0
      ALLOCATE (N_OBSRV_ERR    (BUOYS)); N_OBSRV_ERR = 0
      ALLOCATE (N_STAT         (BUOYS)); N_STAT = 0 

      OPEN(41, FILE = 'scatter.dat', STATUS='UNKNOWN')

!------------------------------------------------------------------------      
! do differences of all integral parameters, all times, all buoys and interpolate in time if needed .... 
!------------------------------------------------------------------------
      
      DO IB = 1, BUOYS

        IF (IFIND(IB) .EQ. 0) THEN

            HS_O(:,IB)               = 0.0
            TM01_O(:,IB)             = 0.0
            TM02_O (:,IB)            = 0.0
            DIFF_HS_SP(:,IB)         = 0.0
            DIFF2_HS_SP(:,IB)        = 0.0
            ABS_DIFF_HS_SP(:,IB)     = 0.0
            DIFF_TM01_SP(:,IB)       = 0.0
            DIFF2_TM01_SP(:,IB)      = 0.0
            ABS_DIFF_TM01_SP (:,IB)  = 0.0
            DIFF_TM02_SP(:,IB)       = 0.0
            DIFF2_TM02_SP(:,IB)      = 0.0
            ABS_DIFF_TM02_SP(:,IB)   = 0.0
            HS_SP_CLEAN(:,IB)        = 0.0
            TM01_SP_CLEAN(:,IB)      = 0.0
            TM02_SP_CLEAN(:,IB)      = 0.0
            HS_O_CLEAN(:,IB)         = 0.0
            TM01_O_CLEAN(:,IB)       = 0.0
            TM02_O_CLEAN(:,IB)       = 0.0

            TP_O(:,IB)               = 0.0
            DIFF_TP_SP(:,IB)         = 0.0
            DIFF2_TP_SP(:,IB)        = 0.0
            ABS_DIFF_TP_SP(:,IB)     = 0.0
            TP_SP_CLEAN(:,IB)        = 0.0
            TP_O_CLEAN(:,IB)         = 0.0

         END IF

        CTMP30 = BNAMES(IB)
        NAMEBUOY = CTMP30(1:5)
        OPEN(40, FILE = NAMEBUOY//'_time_series'//'.dat', STATUS='UNKNOWN')
                    
        COUNTER = 1
        
        DO I = 1, N_DT_O(IB)
 
          !WRITE(1010,'(3I10,3F15.4)') IB, I, N_DT_O(IB), ZEIT_O(I,IB), ZEIT_SP(1), ZEIT_SP(N_DT_SP)
          
          IF (ZEIT_O(I,IB) > ZEIT_SP(N_DT_SP) .OR. ZEIT_O(I,IB) < ZEIT_SP(1)) THEN  ! Wenn die Simulationszeit ausserhalb der Obervationszeit liegt...

            WRITE(1007,'(A50,3I10,3F15.4)') 'Simulationszeit ausserhalb der Obervationszeit', IB, I, N_DT_O(IB), ZEIT_O(I,IB), ZEIT_SP(1), ZEIT_SP(N_DT_SP)
          
            HS_O(I,IB)               = 0.0
            TM01_O(I,IB)             = 0.0
            TM02_O (I,IB)            = 0.0
            DIFF_HS_SP(I,IB)         = 0.0
            DIFF2_HS_SP(I,IB)        = 0.0
            ABS_DIFF_HS_SP(I,IB)     = 0.0
            DIFF_TM01_SP(I,IB)       = 0.0
            DIFF2_TM01_SP(I,IB)      = 0.0
            ABS_DIFF_TM01_SP (I,IB)  = 0.0
            DIFF_TM02_SP(I,IB)       = 0.0
            DIFF2_TM02_SP(I,IB)      = 0.0
            ABS_DIFF_TM02_SP(I,IB)   = 0.0
            HS_SP_CLEAN(I,IB)        = 0.0            
            TM01_SP_CLEAN(I,IB)      = 0.0            
            TM02_SP_CLEAN(I,IB)      = 0.0
            HS_O_CLEAN(I,IB)         = 0.0            
            TM01_O_CLEAN(I,IB)       = 0.0            
            TM02_O_CLEAN(I,IB)       = 0.0           

            TP_O(I,IB)               = 0.0
            DIFF_TP_SP(I,IB)         = 0.0
            DIFF2_TP_SP(I,IB)        = 0.0
            ABS_DIFF_TP_SP(I,IB)     = 0.0
            TP_SP_CLEAN(I,IB)        = 0.0
            TP_O_CLEAN(I,IB)         = 0.0

         
            N_OUT_TIME(IB) = N_OUT_TIME(IB) + 1 

            CYCLE
            
          ELSE
          
            IF (HS_O(I,IB) .LT. SMALL .OR. HS_O(I,IB) .GT. 30.d0) THEN ! Falls Messergebniss falsch ... 
 
              WRITE(1007,'(A50,2I10,2F15.4)') 'Messergebniss falsch ...', IB, I, HS_O(I,IB), HS_O(I,IB) 

              HS_O(I,IB)               = 0.0
              TM01_O(I,IB)             = 0.0
              TM02_O (I,IB)            = 0.0
              DIFF_HS_SP(I,IB)         = 0.0
              DIFF2_HS_SP(I,IB)        = 0.0
              ABS_DIFF_HS_SP(I,IB)     = 0.0
              DIFF_TM01_SP(I,IB)       = 0.0
              DIFF2_TM01_SP(I,IB)      = 0.0
              ABS_DIFF_TM01_SP (I,IB)  = 0.0
              DIFF_TM02_SP(I,IB)       = 0.0
              DIFF2_TM02_SP(I,IB)      = 0.0
              ABS_DIFF_TM02_SP(I,IB)   = 0.0
              HS_SP_CLEAN(I,IB)        = 0.0
              TM01_SP_CLEAN(I,IB)      = 0.0
              TM02_SP_CLEAN(I,IB)      = 0.0
              HS_O_CLEAN(I,IB)         = 0.0
              TM01_O_CLEAN(I,IB)       = 0.0
              TM02_O_CLEAN(I,IB)       = 0.0

              TP_O(I,IB)               = 0.0
              DIFF_TP_SP(I,IB)         = 0.0
              DIFF2_TP_SP(I,IB)        = 0.0
              ABS_DIFF_TP_SP(I,IB)     = 0.0
              TP_SP_CLEAN(I,IB)        = 0.0
              TP_O_CLEAN(I,IB)         = 0.0
            
              N_OBSRV_ERR(IB) = N_OBSRV_ERR(IB) + 1 
 
              CYCLE
              
            ELSE
            
              DO K = 1, N_DT_SP - 1

                IF(ABS(ZEIT_O(I,IB)- ZEIT_SP(K)) .LT. TINY(1.)) THEN  ! Zeit stimmt genau ....

                  WRITE(2003,*) DATUM_O(I,IB), DATUM_SP(K), ZEIT_O(I,IB), ZEIT_SP(K)
                                  
                  HS_SP_CLEAN      (I,IB) =     HS_SP   (K,IB)
                  HS_O_CLEAN       (I,IB) =     HS_O    (I,IB)
                  DIFF_HS_SP       (I,IB) =     HS_SP   (K,IB) -   HS_O (I,IB)
                  DIFF2_HS_SP      (I,IB) =    (HS_SP   (K,IB) -   HS_O (I,IB))**2.
                  ABS_DIFF_HS_SP   (I,IB) = ABS(HS_SP   (K,IB) -   HS_O (I,IB))   
                  
                  TM01_SP_CLEAN    (I,IB) =     TM01_SP (K,IB)
                  TM01_O_CLEAN     (I,IB) =     TM01_O  (I,IB)
                  DIFF_TM01_SP     (I,IB) =     TM01_SP (K,IB) - TM01_O (I,IB)
                  DIFF2_TM01_SP    (I,IB) =    (TM01_SP (K,IB) - TM01_O (I,IB))**2.
                  ABS_DIFF_TM01_SP (I,IB) = ABS(TM01_SP (K,IB) - TM01_O (I,IB))
                  
                  TM02_SP_CLEAN    (I,IB) =     TM02_SP (K,IB)
                  TM02_O_CLEAN     (I,IB) =     TM02_O  (I,IB)
                  DIFF_TM02_SP     (I,IB) =     TM02_SP (K,IB) - TM02_O (I,IB)
                  DIFF2_TM02_SP    (I,IB) =    (TM02_SP (K,IB) - TM02_O (I,IB))**2.
                  ABS_DIFF_TM02_SP (I,IB) = ABS(TM02_SP (K,IB) - TM02_O (I,IB))

                  TP_SP_CLEAN      (I,IB) =     TP_SP   (K,IB)
                  TP_O_CLEAN       (I,IB) =     TP_O    (I,IB)
                  DIFF_TP_SP       (I,IB) =     TP_SP   (K,IB) -   TP_O (I,IB)
                  DIFF2_TP_SP      (I,IB) =    (TP_SP   (K,IB) -   TP_O (I,IB))**2.
                  ABS_DIFF_TP_SP   (I,IB) = ABS(TP_SP   (K,IB) -   TP_O (I,IB))

                  FOUND = .TRUE.
                  COUNTER = COUNTER + 1
                  N_STAT(IB) = N_STAT(IB) + 1 
                  
                  WRITE(40, '(F15.6, 8F15.4)') ZEIT_O(I,IB), HS_SP(K,IB), HS_O(I,IB), &
                                               TM01_SP (K,IB), TM01_O (I,IB), TM02_SP (K,IB), TM02_O (I,IB), TP_SP(K,IB), TP_O(I,IB)
                  WRITE(41, '(F15.6, 8F15.4)') ZEIT_O(I,IB), HS_SP(K,IB), HS_O(I,IB), &
                                               TM01_SP (K,IB), TM01_O (I,IB), TM02_SP (K,IB), TM02_O (I,IB), TP_SP(K,IB), TP_O(I,IB)

                  EXIT 
                  
                ELSE
                
                  FOUND = .FALSE.
                  
                END IF
                              
              END DO
              
              IF (.NOT. FOUND) THEN
              
                DO K = 1, N_DT_SP - 1
                
                  IF(ZEIT_O(I,IB) > ZEIT_SP(K) .AND. ZEIT_O(I,IB) < ZEIT_SP(K+1)) THEN

                    WRITE(2003,*) DATUM_O(I,IB), DATUM_SP(K), ZEIT_O(I,IB), ZEIT_SP(K)

                    INCR_DT  = ZEIT_SP(K+1) - ZEIT_SP(K)
                    INTER_DT = ZEIT_O (I,IB) - ZEIT_SP(K)
                    
                    CALL INTER_S (HS_SP(K,IB),HS_SP(K+1,IB),INCR_DT,INTER_DT,INTER)
                    HS_SP_CLEAN      (I,IB) =     INTER
                    HS_O_CLEAN       (I,IB) =     HS_O    (I,IB)
                    DIFF_HS_SP       (I,IB) =     INTER - HS_O (I,IB)
                    DIFF2_HS_SP      (I,IB) =    (INTER - HS_O (I,IB))**2.
                    ABS_DIFF_HS_SP   (I,IB) = ABS(INTER - HS_O (I,IB)) 
                    
                    CALL INTER_S (TM01_SP (K,IB),TM01_SP (K+1,IB),INCR_DT,INTER_DT,INTER)
                    TM01_SP_CLEAN    (I,IB) =     INTER
                    TM01_O_CLEAN     (I,IB) =     TM01_O  (I,IB)
                    DIFF_TM01_SP     (I,IB) =     INTER - TM01_O (I,IB)
                    DIFF2_TM01_SP    (I,IB) =    (INTER - TM01_O (I,IB))**2.
                    ABS_DIFF_TM01_SP (I,IB) = ABS(INTER - TM01_O (I,IB))  
                    
                    CALL INTER_S (TM02_SP (K,IB),TM02_SP (K+1,IB),INCR_DT,INTER_DT,INTER)
                    TM02_SP_CLEAN    (I,IB) =     INTER
                    TM02_O_CLEAN     (I,IB) =     TM02_O  (I,IB)
                    DIFF_TM02_SP     (I,IB) =     INTER - TM02_O (I,IB)
                    DIFF2_TM02_SP    (I,IB) =    (INTER - TM02_O (I,IB))**2.
                    ABS_DIFF_TM02_SP (I,IB) = ABS(INTER - TM02_O (I,IB)) 

                    CALL INTER_S (TP_SP (K,IB),TP_SP (K+1,IB),INCR_DT,INTER_DT,INTER)
                    TP_SP_CLEAN    (I,IB) =     INTER
                    TP_O_CLEAN     (I,IB) =     TP_O  (I,IB)
                    DIFF_TP_SP     (I,IB) =     INTER - TP_O (I,IB)
                    DIFF2_TP_SP    (I,IB) =    (INTER - TP_O (I,IB))**2.
                    ABS_DIFF_TP_SP (I,IB) = ABS(INTER - TP_O (I,IB))
                    
                    N_STAT(IB) = N_STAT(IB) + 1
                    COUNTER = COUNTER + 1

                    WRITE(40, '(F15.6, 8F15.4)') ZEIT_O(I,IB), HS_SP(K,IB), HS_O(I,IB), &
                                                  TM01_SP (K,IB), TM01_O (I,IB), TM02_SP (K,IB), TM02_O (I,IB), TP_SP(K,IB), TP_O(I,IB)
                    WRITE(41, '(F15.6, 8F15.4)') ZEIT_O(I,IB), HS_SP(K,IB), HS_O(I,IB), &
                                                  TM01_SP (K,IB), TM01_O (I,IB), TM02_SP (K,IB), TM02_O (I,IB), TP_SP(K,IB), TP_O(I,IB)
 
                    EXIT

                  END IF
                END DO              
              END IF
            END IF ! HS_O(I,IB) .LT. SMALL .OR. HS_O(I,IB) .GT. 30.d0
          END IF ! ZEIT_O(I,IB) > ZEIT_SP(N_DT_SP) .OR. ZEIT_O(I,IB) < ZEIT_SP(1)
        END DO ! I
        CLOSE(40)
      END DO ! IB

!------------------------------------------------------------------------      
! allocate arrays for mean values ... 
!------------------------------------------------------------------------
       
      ALLOCATE (MEAN_HS_O   (BUOYS))
      ALLOCATE (BIAS_HS_SP  (BUOYS))
      ALLOCATE (MAE_HS_SP   (BUOYS))
      ALLOCATE (RMS_HS_SP   (BUOYS))
      ALLOCATE (SCI_HS_SP   (BUOYS))
      ALLOCATE (MEAN_TM01_O (BUOYS))
      ALLOCATE (BIAS_TM01_SP(BUOYS))
      ALLOCATE (MAE_TM01_SP (BUOYS))
      ALLOCATE (RMS_TM01_SP (BUOYS))
      ALLOCATE (SCI_TM01_SP (BUOYS))
      ALLOCATE (MEAN_TM02_O (BUOYS))
      ALLOCATE (BIAS_TM02_SP(BUOYS))
      ALLOCATE (MAE_TM02_SP (BUOYS))
      ALLOCATE (RMS_TM02_SP (BUOYS))    
      ALLOCATE (SCI_TM02_SP (BUOYS))       
      ALLOCATE (MEAN_HS_SP  (BUOYS))
      ALLOCATE (MEAN_TM01_SP(BUOYS))
      ALLOCATE (MEAN_TM02_SP(BUOYS))
      ALLOCATE (KORR_HS     (BUOYS))
      ALLOCATE (KORR_TM01   (BUOYS))
      ALLOCATE (KORR_TM02   (BUOYS))  

      ALLOCATE (MEAN_TP_O   (BUOYS))
      ALLOCATE (BIAS_TP_SP  (BUOYS))
      ALLOCATE (MAE_TP_SP   (BUOYS))
      ALLOCATE (RMS_TP_SP   (BUOYS))
      ALLOCATE (SCI_TP_SP   (BUOYS))
      ALLOCATE (MEAN_TP_SP  (BUOYS))
      ALLOCATE (KORR_TP     (BUOYS))

      
      MEAN_HS_O   (:) = 0.0
      BIAS_HS_SP  (:) = 0.0
      MAE_HS_SP   (:) = 0.0
      RMS_HS_SP   (:) = 0.0
      SCI_HS_SP   (:) = 0.0
      MEAN_TM01_O (:) = 0.0
      BIAS_TM01_SP(:) = 0.0
      MAE_TM01_SP (:) = 0.0
      RMS_TM01_SP (:) = 0.0
      SCI_TM01_SP (:) = 0.0
      MEAN_TM02_O (:) = 0.0
      BIAS_TM02_SP(:) = 0.0
      MAE_TM02_SP (:) = 0.0
      RMS_TM02_SP (:) = 0.0
      SCI_TM02_SP (:) = 0.0   
      MEAN_HS_SP  (:) = 0.0
      MEAN_TM01_SP(:) = 0.0
      MEAN_TM02_SP(:) = 0.0
      KORR_HS     (:) = 0.0
      KORR_TM01   (:) = 0.0
      KORR_TM02   (:) = 0.0    

      MEAN_TP_O   (:) = 0.0
      BIAS_TP_SP  (:) = 0.0
      MAE_TP_SP   (:) = 0.0
      RMS_TP_SP   (:) = 0.0
      SCI_TP_SP   (:) = 0.0
      MEAN_TP_SP  (:) = 0.0
      KORR_TP     (:) = 0.0

      !WRITE(*,'(A10,A10,5A16)') 'BUOY NO', 'BNAMES', 'IFIND', 'N_STAT', 'N_OUT', 'N_OBSRV_ERR'  
      !DO I = 1, BUOYS
      !  WRITE (*,'(I10,A10,5I15)') I, TRIM(BNAMES(I)), IFIND(I), N_STAT(I), N_OUT_TIME(I), N_OBSRV_ERR(I)
      !END DO  

!------------------------------------------------------------------------      
! read the observed wave spectra and estimate spectral moments ...
!------------------------------------------------------------------------

      DO I = 1, BUOYS
      
        IF (IFIND(I) .EQ. 0 .OR. N_STAT(I) .EQ. 0) THEN

          WRITE(1007,*) 'BUOY NOT FOUND NOR ANY VALUABLE DATA', IFIND(I), N_STAT(I)

          MEAN_HS_O (I)   = 0.
          MEAN_HS_SP(I)   = 0.
          BIAS_HS_SP(I)   = 0.
          MAE_HS_SP (I)   = 0.
          RMS_HS_SP (I)   = 0.
          SCI_HS_SP (I)   = 0.
        
          MEAN_TM01_O (I) = 0.
          MEAN_TM01_SP(I) = 0.
          BIAS_TM01_SP(I) = 0.
          MAE_TM01_SP (I) = 0.
          RMS_TM01_SP (I) = 0.
          SCI_TM01_SP (I) = 0.
          
          MEAN_TM02_O (I) = 0.
          MEAN_TM02_SP(I) = 0.
          BIAS_TM02_SP(I) = 0.
          MAE_TM02_SP (I) = 0.
          RMS_TM02_SP (I) = 0.
          SCI_TM02_SP (I) = 0.

          MEAN_TP_O (I) = 0.
          MEAN_TP_SP(I) = 0.
          BIAS_TP_SP(I) = 0.
          MAE_TP_SP (I) = 0.
          RMS_TP_SP (I) = 0.
          SCI_TP_SP (I) = 0.

        ELSE
 
          MEAN_HS_O (I)   =  SUM(HS_O_CLEAN(:,I))      / (N_STAT(I))
          MEAN_HS_SP(I)   =  SUM(HS_SP_CLEAN(:,I))     / (N_STAT(I))
          BIAS_HS_SP(I)   =  SUM(DIFF_HS_SP    (:,I))  / (N_STAT(I))
          MAE_HS_SP (I)   =  SUM(ABS_DIFF_HS_SP(:,I))  / (N_STAT(I))
          RMS_HS_SP (I)   =  (SUM((DIFF_HS_SP(:,I)-BIAS_HS_SP(I))**2.)/N_STAT(I))**0.5
          SCI_HS_SP (I)   =  RMS_HS_SP (I) / MEAN_HS_O (I)
        
          MEAN_TM01_O (I) =  SUM(TM01_O_CLEAN(:,I))     / (N_STAT(I))
          MEAN_TM01_SP(I) =  SUM(TM01_SP_CLEAN(:,I))    / (N_STAT(I))
          BIAS_TM01_SP(I) =  SUM(DIFF_TM01_SP    (:,I)) / (N_STAT(I))
          MAE_TM01_SP (I) =  SUM(ABS_DIFF_TM01_SP(:,I)) / (N_STAT(I))
          RMS_TM01_SP (I)   =  (SUM((DIFF_TM01_SP(:,I)-BIAS_TM01_SP(I))**2.)/N_STAT(I))**0.5
          SCI_TM01_SP (I) =  RMS_TM01_SP (I) / MEAN_TM01_O (I)        
          
          MEAN_TM02_O (I) =  SUM(TM02_O_CLEAN(:,I))     / (N_STAT(I))
          MEAN_TM02_SP(I) =  SUM(TM02_SP_CLEAN(:,I))    / (N_STAT(I))
          BIAS_TM02_SP(I) =  SUM(DIFF_TM02_SP    (:,I)) / (N_STAT(I))
          MAE_TM02_SP (I) =  SUM(ABS_DIFF_TM02_SP(:,I)) / (N_STAT(I))
          RMS_TM02_SP (I) =  (SUM((DIFF_TM02_SP(:,I)-BIAS_TM02_SP(I))**2.)/N_STAT(I))**0.5
          SCI_TM02_SP (I) =  RMS_TM02_SP (I) / MEAN_TM02_O (I)

          MEAN_TP_O (I) =  SUM(TP_O_CLEAN(:,I))     / (N_STAT(I))
          MEAN_TP_SP(I) =  SUM(TP_SP_CLEAN(:,I))    / (N_STAT(I))
          BIAS_TP_SP(I) =  SUM(DIFF_TP_SP    (:,I)) / (N_STAT(I))
          MAE_TP_SP (I) =  SUM(ABS_DIFF_TP_SP(:,I)) / (N_STAT(I))
          RMS_TP_SP (I) =  (SUM((DIFF_TP_SP(:,I)-BIAS_TP_SP(I))**2.)/N_STAT(I))**0.5
          SCI_TP_SP (I) =  RMS_TP_SP (I) / MEAN_TP_O (I)


        END IF
        
      END DO
      
      MEAN_O_HS_ALL    = SUM(HS_O_CLEAN(:,:)) /  SUM(N_STAT(:))
      MEAN_SP_HS_ALL   = SUM(HS_SP_CLEAN(:,:)) /  SUM(N_STAT(:))
      MEAN_O_TM01_ALL  = SUM(TM01_O_CLEAN(:,:)) /  SUM(N_STAT(:))
      MEAN_SP_TM01_ALL = SUM(TM01_SP_CLEAN(:,:)) /  SUM(N_STAT(:))
      MEAN_O_TM02_ALL  = SUM(TM02_O_CLEAN(:,:)) /  SUM(N_STAT(:))
      MEAN_SP_TM02_ALL = SUM(TM02_SP_CLEAN(:,:)) /  SUM(N_STAT(:))
      MEAN_O_TP_ALL    = SUM(TP_O_CLEAN(:,:)) /  SUM(N_STAT(:))
      MEAN_SP_TP_ALL   = SUM(TP_SP_CLEAN(:,:)) /  SUM(N_STAT(:))
      
      BIAS_HS_ALL   = SUM(DIFF_HS_SP(:,:)) /  SUM(N_STAT(:))
      BIAS_TM01_ALL = SUM(DIFF_TM01_SP(:,:)) /  SUM(N_STAT(:))
      BIAS_TM02_ALL = SUM(DIFF_TM02_SP(:,:)) /  SUM(N_STAT(:))
      BIAS_TP_ALL   = SUM(DIFF_TP_SP(:,:)) /  SUM(N_STAT(:))
      
      RMS_HS_ALL   = (SUM(DIFF2_HS_SP(:,:)) /  SUM(N_STAT(:)))**0.5
      RMS_TM01_ALL = (SUM(DIFF2_TM01_SP(:,:)) /  SUM(N_STAT(:)))**0.5
      RMS_TM02_ALL = (SUM(DIFF2_TM02_SP(:,:)) /  SUM(N_STAT(:)))**0.5
      RMS_TP_ALL = (SUM(DIFF2_TP_SP(:,:)) /  SUM(N_STAT(:)))**0.5
      
      SCI_HS_ALL   = RMS_HS_ALL   / MEAN_O_HS_ALL
      SCI_TM01_ALL = RMS_TM01_ALL / MEAN_O_TM01_ALL
      SCI_TM02_ALL = RMS_TM02_ALL / MEAN_O_TM02_ALL  
      SCI_TP_ALL = RMS_TP_ALL / MEAN_O_TP_ALL
      
      MAE_HS_ALL   = SUM(ABS_DIFF_HS_SP  (:,:)) /  SUM(N_STAT(:))
      MAE_TM01_ALL = SUM(ABS_DIFF_TM01_SP(:,:)) /  SUM(N_STAT(:))
      MAE_TM02_ALL = SUM(ABS_DIFF_TM02_SP(:,:)) /  SUM(N_STAT(:)) 
      MAE_TP_ALL = SUM(ABS_DIFF_TP_SP(:,:)) /  SUM(N_STAT(:))

!------------------------------------------------------------------------      
! calculate the differenceses for the correlation coefficients ...
!------------------------------------------------------------------------
 
      DIFF_HS_SP = 0.
      DIFF2_HS_SP = 0.
      DIFF3_HS_SP = 0.
      DIFF4_HS_SP = 0.

      DO IB = 1, BUOYS

        IF (IFIND(IB) .GT. 0 .AND. N_STAT(IB) .GT. 0 ) THEN

          DO I = 1, N_DT_O(IB)

            IF (HS_O_CLEAN    (I,IB) .GT. SMALL .AND. HS_O_CLEAN    (I,IB) .LT. 30.)THEN

              DIFF_HS_SP    (I,IB) =  HS_O_CLEAN    (I,IB) - MEAN_HS_O    (IB)
              DIFF2_HS_SP   (I,IB) = (HS_O_CLEAN    (I,IB) - MEAN_HS_O    (IB))**2.0
              DIFF3_HS_SP   (I,IB) =  HS_SP_CLEAN   (I,IB) - MEAN_HS_SP   (IB)
              DIFF4_HS_SP   (I,IB) = (HS_SP_CLEAN   (I,IB) - MEAN_HS_SP   (IB))**2.0
            
              DIFF_TM01_SP  (I,IB) =  TM01_O_CLEAN  (I,IB) - MEAN_TM01_O  (IB)
              DIFF2_TM01_SP (I,IB) = (TM01_O_CLEAN  (I,IB) - MEAN_TM01_O  (IB))**2.0
              DIFF3_TM01_SP (I,IB) =  TM01_SP_CLEAN (I,IB) - MEAN_TM01_SP (IB)
              DIFF4_TM01_SP (I,IB) = (TM01_SP_CLEAN (I,IB) - MEAN_TM01_SP (IB))**2.0  
            
              DIFF_TM02_SP  (I,IB) =  TM02_O_CLEAN  (I,IB) - MEAN_TM02_O  (IB)
              DIFF2_TM02_SP (I,IB) = (TM02_O_CLEAN  (I,IB) - MEAN_TM02_O  (IB))**2.0
              DIFF3_TM02_SP (I,IB) =  TM02_SP_CLEAN (I,IB) - MEAN_TM02_SP (IB)
              DIFF4_TM02_SP (I,IB) = (TM02_SP_CLEAN (I,IB) - MEAN_TM02_SP (IB))**2.0

              DIFF_TP_SP  (I,IB) =  TP_O_CLEAN  (I,IB) - MEAN_TP_O  (IB)
              DIFF2_TP_SP (I,IB) = (TP_O_CLEAN  (I,IB) - MEAN_TP_O  (IB))**2.0
              DIFF3_TP_SP (I,IB) =  TP_SP_CLEAN (I,IB) - MEAN_TP_SP (IB)
              DIFF4_TP_SP (I,IB) = (TP_SP_CLEAN (I,IB) - MEAN_TP_SP (IB))**2.0


            ELSE

              DIFF_HS_SP    (I,IB) = 0.
              DIFF2_HS_SP   (I,IB) = 0.
              DIFF3_HS_SP   (I,IB) = 0.
              DIFF4_HS_SP   (I,IB) = 0.
  
              DIFF_TM01_SP  (I,IB) = 0.
              DIFF2_TM01_SP (I,IB) = 0.
              DIFF3_TM01_SP (I,IB) = 0.
              DIFF4_TM01_SP (I,IB) = 0.

              DIFF_TM02_SP  (I,IB) = 0.
              DIFF2_TM02_SP (I,IB) = 0.
              DIFF3_TM02_SP (I,IB) = 0.
              DIFF4_TM02_SP (I,IB) = 0.

              DIFF_TP_SP  (I,IB) = 0.
              DIFF2_TP_SP (I,IB) = 0.
              DIFF3_TP_SP (I,IB) = 0.
              DIFF4_TP_SP (I,IB) = 0.

            END IF

          END DO

        ELSE

            DIFF_HS_SP    (:,IB) = 0. 
            DIFF2_HS_SP   (:,IB) = 0.
            DIFF3_HS_SP   (:,IB) = 0.
            DIFF4_HS_SP   (:,IB) = 0.

            DIFF_TM01_SP  (:,IB) = 0.
            DIFF2_TM01_SP (:,IB) = 0.
            DIFF3_TM01_SP (:,IB) = 0.
            DIFF4_TM01_SP (:,IB) = 0.

            DIFF_TM02_SP  (:,IB) = 0. 
            DIFF2_TM02_SP (:,IB) = 0.
            DIFF3_TM02_SP (:,IB) = 0.
            DIFF4_TM02_SP (:,IB) = 0. 

            DIFF_TP_SP  (:,IB) = 0.
            DIFF2_TP_SP (:,IB) = 0.
            DIFF3_TP_SP (:,IB) = 0.
            DIFF4_TP_SP (:,IB) = 0.


         END IF
        
      END DO  

!------------------------------------------------------------------------      
! calculate correlation coefficients ...
!------------------------------------------------------------------------

      DO I = 1, BUOYS
        IF (IFIND(I) .EQ. 0 .OR. N_STAT(I) .EQ. 0) THEN
           KORR_HS(I) = 0.
           KORR_TM01(I) = 0.
           KORR_TM02(I) = 0.
           KORR_TP(I) = 0.
        ELSE
           KORR_HS(I)   = (SUM(DIFF_HS_SP(:,I)   * DIFF3_HS_SP(:,I))   / &
                      ((SUM(DIFF2_HS_SP(:,I))   * SUM(DIFF4_HS_SP(:,I)))   **0.5))!**2.
           KORR_TM01(I) = (SUM(DIFF_TM01_SP(:,I) * DIFF3_TM01_SP(:,I)) / & 
                      ((SUM(DIFF2_TM01_SP(:,I)) * SUM(DIFF4_TM01_SP(:,I))) **0.5))!**2.
           KORR_TM02(I) = (SUM(DIFF_TM02_SP(:,I) * DIFF3_TM02_SP(:,I)) / & 
                      ((SUM(DIFF2_TM02_SP(:,I)) * SUM(DIFF4_TM02_SP(:,I))) **0.5))!**2.
           KORR_TP(I) = (SUM(DIFF_TP_SP(:,I) * DIFF3_TP_SP(:,I)) / &
                      ((SUM(DIFF2_TP_SP(:,I)) * SUM(DIFF4_TP_SP(:,I))) **0.5))!**2.

        END IF
      END DO
    
      IF (SUM(ABS(DIFF_HS_SP(:,:))) .GT. SMALL) THEN 
        KORR_HS_ALL   =  (SUM(DIFF_HS_SP(:,:)   * DIFF3_HS_SP(:,:))   / &
                        ((SUM(DIFF2_HS_SP(:,:))    * SUM(DIFF4_HS_SP(:,:)))    **0.5))!**2.
        KORR_TM01_ALL =  (SUM(DIFF_TM01_SP(:,:) * DIFF3_TM01_SP(:,:)) / &
                        ((SUM(DIFF2_TM01_SP(:,:))  * SUM(DIFF4_TM01_SP(:,:)))  **0.5))!**2.
        KORR_TM02_ALL =  (SUM(DIFF_TM02_SP(:,:) * DIFF3_TM02_SP(:,:)) / & 
                        ((SUM(DIFF2_TM02_SP(:,:))  * SUM(DIFF4_TM02_SP(:,:)))  **0.5))!**2.
        KORR_TM02_ALL =  (SUM(DIFF_TP_SP(:,:) * DIFF3_TP_SP(:,:)) / &
                        ((SUM(DIFF2_TP_SP(:,:))  * SUM(DIFF4_TP_SP(:,:)))  **0.5))!**2.

      ELSE
        KORR_HS_ALL   =  0.!
        KORR_TM01_ALL =  0.!
        KORR_TM02_ALL =  0.!
        KORR_TP_ALL =  0.!
      END IF
      
      WRITE(TMPCHAR,999) BUOYS 
     
      FMT_NAM_HEADER = '('//'A12,'//TMPCHAR//'A8'//')' 
      FMT_ERG_HEADER = '('//'A12,'//TMPCHAR//'A8'//')'
      FMT_ERG_RESULT = '('//'A12,'//TMPCHAR//'F8.3'//')'
        
!      WRITE(*,*) FMT_ERG_HEADER
!      WRITE(*,*) FMT_ERG_RESULT

!------------------------------------------------------------------------      
! write statistics of the spectral moments ... 
!------------------------------------------------------------------------
      
      OPEN(10,  FILE = 'statistik.dat', STATUS='UNKNOWN')

      WRITE (10,*) '-----------------------------------------'
      WRITE (10,FMT_NAM_HEADER) 'BUOYS ', (BNAMES(I),I = 1, BUOYS)
      WRITE (10,*) '-----------------HS----------------------'
      WRITE (10,FMT_ERG_RESULT) 'MEAN_O', (MEAN_HS_O (I)  ,I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'MEAN_S', (MEAN_HS_SP(I)  ,I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'BIAS'  , (BIAS_HS_SP(I)  ,I = 1, BUOYS) 
      WRITE (10,FMT_ERG_RESULT) 'MEA'  ,  (MAE_HS_SP (I)  ,I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'RMS'  ,  (RMS_HS_SP (I)  ,I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'SCI'  ,  (SCI_HS_SP (I)  ,I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'KORR' ,  (KORR_HS   (I)  ,I = 1, BUOYS)
      WRITE (10,*) '-----------------TM01--------------------' 
      WRITE (10,FMT_ERG_RESULT) 'MEAN_O', (MEAN_TM01_O (I),I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'MEAN_S', (MEAN_TM01_SP(I),I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'BIAS'  , (BIAS_TM01_SP(I),I = 1, BUOYS) 
      WRITE (10,FMT_ERG_RESULT) 'MEA'  ,  (MAE_TM01_SP (I),I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'RMS'  ,  (RMS_TM01_SP (I),I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'SCI'  ,  (SCI_TM01_SP (I),I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'KORR' ,  (KORR_TM01   (I),I = 1, BUOYS)
      WRITE (10,*) '-----------------TM02--------------------'
      WRITE (10,FMT_ERG_RESULT) 'MEAN_O', (MEAN_TM02_O (I),I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'MEAN_S', (MEAN_TM02_SP(I),I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'BIAS'  , (BIAS_TM02_SP(I),I = 1, BUOYS) 
      WRITE (10,FMT_ERG_RESULT) 'MEA'  ,  (MAE_TM02_SP (I),I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'RMS'  ,  (RMS_TM02_SP (I),I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'SCI'  ,  (SCI_TM02_SP (I),I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'KORR' ,  (KORR_TM02   (I),I = 1, BUOYS)
      WRITE (10,*) '-----------------TP--------------------'
      WRITE (10,FMT_ERG_RESULT) 'MEAN_O', (MEAN_TP_O (I),I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'MEAN_S', (MEAN_TP_SP(I),I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'BIAS'  , (BIAS_TP_SP(I),I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'MEA'  ,  (MAE_TP_SP (I),I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'RMS'  ,  (RMS_TP_SP (I),I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'SCI'  ,  (SCI_TP_SP (I),I = 1, BUOYS)
      WRITE (10,FMT_ERG_RESULT) 'KORR' ,  (KORR_TP   (I),I = 1, BUOYS)


      WRITE (10,*) 'TOTAL SCORE'
      
      FMT_ERG_ALL_STAT = '(A8, F10.3, A8, F10.3)'
      
      WRITE (10,*) '------------------HS------------------'
      WRITE (10,FMT_ERG_ALL_STAT) 'MEAN_O', MEAN_O_HS_ALL,   'MEAN_SP', MEAN_SP_HS_ALL 
      WRITE (10,FMT_ERG_ALL_STAT) 'BIAS'  , BIAS_HS_ALL,     'RMS'    , RMS_HS_ALL 
      WRITE (10,FMT_ERG_ALL_STAT) 'SCI'   , SCI_HS_ALL,      'KORR'   , KORR_HS_ALL 
      WRITE (10,*) '------------------TM01------------------'
      WRITE (10,FMT_ERG_ALL_STAT) 'MEAN_O', MEAN_O_TM01_ALL, 'MEAN_SP', MEAN_SP_TM01_ALL 
      WRITE (10,FMT_ERG_ALL_STAT) 'BIAS'  , BIAS_TM01_ALL,   'RMS'    , RMS_TM01_ALL 
      WRITE (10,FMT_ERG_ALL_STAT) 'SCI'   , SCI_TM01_ALL,    'KORR'   , KORR_TM01_ALL 
      WRITE (10,*) '------------------TM02------------------' 
      WRITE (10,FMT_ERG_ALL_STAT) 'MEAN_O', MEAN_O_TM02_ALL, 'MEAN_SP', MEAN_SP_TM02_ALL 
      WRITE (10,FMT_ERG_ALL_STAT) 'BIAS'  , BIAS_TM02_ALL,   'RMS'    , RMS_TM02_ALL 
      WRITE (10,FMT_ERG_ALL_STAT) 'SCI'   , SCI_TM02_ALL,    'KORR'   , KORR_TM02_ALL  
      WRITE (10,*) '------------------TP------------------'
      WRITE (10,FMT_ERG_ALL_STAT) 'MEAN_O', MEAN_O_TP_ALL, 'MEAN_SP', MEAN_SP_TP_ALL
      WRITE (10,FMT_ERG_ALL_STAT) 'BIAS'  , BIAS_TP_ALL,   'RMS'    , RMS_TP_ALL
      WRITE (10,FMT_ERG_ALL_STAT) 'SCI'   , SCI_TP_ALL,    'KORR'   , KORR_TP_ALL
      
      CLOSE (10)

!------------------------------------------------------------------------      
! start with spectral statistics ...  
!------------------------------------------------------------------------

      ALLOCATE (E_B(N_DT_SP, MSC))

      ALLOCATE (DIFF1_ESP(N_DT_O_MAX, MSC))
      ALLOCATE (DIFF2_ESP(N_DT_O_MAX, MSC))
      ALLOCATE (DIFF1_K(N_DT_O_MAX, MSC))
      ALLOCATE (DIFF2_K(N_DT_O_MAX, MSC))
      ALLOCATE (DIFF3_K(N_DT_O_MAX, MSC))
      ALLOCATE (DIFF4_K(N_DT_O_MAX, MSC))
      ALLOCATE (MEAN_BSP(BUOYS,MSC))
      ALLOCATE (MEAN_OSP(BUOYS,MSC))
      ALLOCATE (BIAS_SP(BUOYS,MSC))
      ALLOCATE (KORR_SP(BUOYS,MSC))
      ALLOCATE (N_DT_ERR(BUOYS, N_DT_O_MAX))
      
      DIFF1_ESP(:,:) = 0.0
      DIFF2_ESP(:,:) = 0.0
      DIFF1_K  (:,:) = 0.0
      DIFF2_K  (:,:) = 0.0
      DIFF3_K  (:,:) = 0.0
      DIFF4_K  (:,:) = 0.0
      MEAN_BSP (:,:) = 0.0
      MEAN_OSP (:,:) = 0.0
      BIAS_SP  (:,:) = 0.0
      KORR_SP  (:,:) = 0.0
      N_DT_ERR(:,:) = 0

      OPEN(101, FILE = 'specstat.dat', STATUS='UNKNOWN')
      
      DO IB = 2, BUOYS
        IF (MAXVAL(FREQ_O(:,IB)) .LT. MAXVAL(FREQ_O(:,IB-1))) MINMAXFREQ_O = MAXVAL(FREQ_O(:,IB))
      END DO      
      
      DO IB = 1, BUOYS
        
        IF (FILETYPE(IB) .EQ. 2) CYCLE

        IF ((IFIND(IB) .EQ. 0).or.(N_STAT(IB) .EQ. 0)) CYCLE

        CTMP30 = BNAMES(IB)
        NAMEBUOY = CTMP30(1:5)

        CALL TEST_FILE_EXIST_DIE("Missing station file : ", TRIM(STATIONNAMES(IB)))
        OPEN(4, FILE = TRIM(STATIONNAMES(IB)), STATUS='OLD', FORM='UNFORMATTED')
        READ(4) MSC
        READ(4) MDC
        READ(4) SPSIG
        READ(4) SPDIR
        READ(4) IFIND(IB)
        DO I = 1, N_DT_SP
          READ  (4) DATUM_SP(I) 
          !READ(4) WINDXY(1), WINDXY(2)
          READ  (4) DEPLOC
          READ  (4) CURTXYLOC
          CALL CT2MJD(DATUM_SP(I), XMJD)
          ZEIT_SP(I) = XMJD 
          READ(4) ACLOC
          READ(4) ACLOC2D
          !WRITE(1008,*) SUM(ACLOC), SUM(ACLOC2D)
          DO IS = 1, MSC
            E_B(I,IS) = SUM(ACLOC2D(IS,:)) * DDIR * PI2 !* SPSIG(IS) 
          END DO
          ETOT = 0 
!          DO IS = 1, MSC - 1
!            ETOT = ETOT + E_B(I,IS) * (SPSIG(IS+1)-SPSIG(IS))/PI2
!          END DO 
!          WRITE(*,*) 'HS ---------->', 4*SQRT(ETOT)
        END DO
        CLOSE(4) 
       
        !WRITE(*,*) STATIONNAMES_O(IB), N_DT_O(IB)

        N_DT_ERR = 0
        CALL TEST_FILE_EXIST_DIE("Missing obse file : ", TRIM(STATIONNAMES_O(IB)))
        OPEN(5, FILE = TRIM(STATIONNAMES_O(IB)), STATUS='OLD')  
        READ (5, '(A500)') DUMP 
        DO IT = 1, N_DT_O(IB)
          ALLOCATE (TMPE(MSC_O(IB)))
          IF(dateFMT == 1) THEN
!           YYYY MM DD hh
            READ  (5, '(A4,A,A2,A,A2,A,A2'//','//FMT_O(IB)//')') YYYY,DUMP3,MM,DUMP3,DD,DUMP3,hh, TMPE(:)
            !WRITE(2222, '(A4,A,A2,A,A2,A,A2)') YYYY,DUMP3,MM,DUMP3,DD,DUMP3,hh
            DATUM_O(IT,IB) = YYYY//MM//DD//'.'//HH//'0000'
          ELSE IF(dateFMT == 2) THEN
!             YYYY MM DD hh MINUTE
            READ  (5, '(A4,A,A2,A,A2,A,A2,A,A2'//','//FMT_O(IB)//')') YYYY,DUMP3,MM,DUMP3,DD,DUMP3,hh,DUMP3,MINUTE, TMPE(:)
            DATUM_O(IT,IB) = YYYY//MM//DD//'.'//HH//MINUTE//'00'
          END IF
          DO IS_O = 1, MSC_O(IB)
            IF (TMPE(IS_O) .GT. 990. .OR. TMPE(IS_O) .LT. 0.) THEN
              N_DT_ERR(IB,IT) = 1 
              EXIT
            ELSE
              N_DT_ERR(IB,IT) = 0 
            END IF 
            IF (FREQ_O(IS_O,IB) .LE. CUT_OFF) THEN
              E_OS(IT,IS_O) = TMPE(IS_O)
            ELSE
              E_OS(IT,IS_O) = 0. 
            END IF
          END DO
!         WRITE (*,*) MSC_O(IB), FMT_O(IB)
          DEALLOCATE(TMPE)
!         WRITE (*,*) 'ZEIT', I, N_DT_O(IB)
         !WRITE  (*, '(A14'//','//FMT_O//')') DATUM_O(I,IB), E_OS(I,:)
          CALL CT2MJD(DATUM_O(IT,IB), XMJD)
          ZEIT_O(IT,IB) = XMJD
        END DO
        CLOSE(5)        
!        
!       Interpolation der Observationsergebnisse auf Frequenzen der Berechnungen (Bojenmessungen meißt dichter im Frequenzraum)               
!       Volle Schleife über die Frequenzen der Berechnungen 
!
        DO IT_O = 1, N_DT_O(IB)
          DO IS = 1, MSC
            IF (FREQ_SP(IS) .LE. FREQ_O(MSC_O(IB),IB)) THEN
              DO IS_O = 1, MSC_O(IB) - 1
                IF ( FREQ_SP(IS) .GE. FREQ_O(IS_O,IB) .AND. FREQ_SP(IS) .LT. FREQ_O(IS_O+1,IB) ) THEN
                  DX      =  FREQ_O(IS_O+1,IB) - FREQ_O(IS_O,IB)
                  DIFF_DX =  FREQ_SP(IS)       - FREQ_O(IS_O,IB) 
                  CALL INTER_S ( E_OS(IT_O,IS_O) , E_OS(IT_O,IS_O+1) ,DX, DIFF_DX, YINTER)
                  E_OSF(IT_O,IS) = YINTER
                  IF (YINTER .LT. 0.) THEN
                    WRITE(*,'(9F10.5)') E_OS(IT_O,IS_O), E_OS(IT_O,IS_O+1), FREQ_SP(IS), FREQ_SP(IS+1), FREQ_O(IS_O,IB), FREQ_O(IS_O+1,IB), DX, DIFF_DX, YINTER 
                    STOP 'ERROR IN FREQ. INTERPOLATION'
                  END IF
                END IF
              END DO
            ELSE
              E_OSF(IT_O,IS) = 0.0
            END IF
          END DO
        END DO 
!
!       Interpolation der Berechnungsergebnisse auf die Zeiten der Observationen (Berechnungsergebnisse meißt dichter im Zeitraum )        
!       Volle Schleife über die Zeiten der Observationen      
!
        DO IT_O = 1, N_DT_O(IB)
          !WRITE(*,'(4F15.4)') ZEIT_O(IT,IB), ZEIT_SP(N_DT_SP), ZEIT_SP(1)
          IF (ZEIT_O(IT_O,IB) > ZEIT_SP(1) .AND. ZEIT_O(IT_O,IB) < ZEIT_SP(N_DT_SP)) THEN ! Check if the observation time is within the simulation period
            !WRITE(*,'(4F15.4)') ZEIT_O(IT_O,IB), ZEIT_SP(1), ZEIT_SP(N_DT_SP)
            DO IT = 1, N_DT_SP - 1  
              IF (ZEIT_O(IT_O,IB) .GE. ZEIT_SP(IT) .AND. ZEIT_O(IT_O,IB) .LT. ZEIT_SP(IT+1)) THEN
                DO IS = 1, MSC
                  IF (FREQ_SP(IS) .LE. FREQ_O(MSC_O(IB),IB)) THEN
                    DT      =  ZEIT_SP(IT+1) - ZEIT_SP(IT)
                    DIFF_DT =  ZEIT_O(IT_O,IB) - ZEIT_SP(IT)  
                    CALL INTER_S (E_B(IT,IS), E_B(IT+1,IS), DT, DIFF_DT, YINTER)
                    E_BT(IT_O,IS) = YINTER
                    !WRITE(1009,*) IS, FREQ_SP(IS), E_BT(I,IS)
                  ELSE
                    E_BT(IT_O,IS) = 0.0
                  END IF
                  IF (YINTER .LT. 0.0) THEN
                    WRITE (*,*) ZEIT_SP(IT_O), ZEIT_SP(IT_O+1), ZEIT_O(IT,IB), ZEIT_O(IT+1,IB)
                    WRITE (*,'(6F10.5)') FREQ_SP(IS), E_B(IT,IS), E_B(IT+1,IS), YINTER, DT, DIFF_DT
                    STOP 'ERROR IN TIME INTERPOLATION'
                  END IF 
                END DO
              END IF ! Time Window
            END DO
          ELSE
            E_BT(IT_O,:) = 0.
          END IF
        END DO        
!        
!       Berechnen der Differenzen zwischen berechneten und gemessenen Energiedifferenzen
!    
        N_STAT(IB) = 0
        OPEN(48, FILE = TRIM(NAMEBUOY)//'_diffspec'//'.dat', STATUS='UNKNOWN')
        OPEN(49, FILE = TRIM(NAMEBUOY)//'_spec'//'.dat', STATUS='UNKNOWN')
!        OPEN(50, FILE = TRIM(NAMEBUOY)//'_observation'//'.dat', STATUS='UNKNOWN')
        WRITE(CHTMP,999) MSC       
        !WRITE(48,'('//CHTMP//'F9.2'//')') (FREQ_SP(I), I = 1, MSC) 
        DO IT_O = 1, N_DT_O(IB) 
          IF (ZEIT_O(IT_O,IB) > ZEIT_SP(1) .AND. ZEIT_O(IT_O,IB) < ZEIT_SP(N_DT_SP)) THEN
            !write(*,*) IT_O, N_DT_O(IB), ZEIT_O(IT_O,IB), ZEIT_SP(1), ZEIT_O(IT_O,IB), ZEIT_SP(N_DT_SP)
            IF (N_DT_ERR(IB,IT_O) .EQ. 0) THEN
              ETMPS = 0.
              ETMPO = 0.
              N_STAT(IB) = N_STAT(IB) + 1 
              !write(*,*) IB
              DO IS = 1, MSC
                !write(*,*) FREQ_SP(IS), FREQ_O(MSC_O(IB),IB)
                IF((FREQ_SP(IS) .LE. FREQ_O(MSC_O(IB),IB))) THEN
                  DIFF1_ESP(IT_O,IS) =  E_BT(IT_O,IS) - E_OSF(IT_O,IS)
                  DIFF2_ESP(IT_O,IS) =  DIFF1_ESP(IT_O,IS)**2. 
                  WRITE(48,'(I10,6F20.8)') IT_O, ZEIT_O(IT_O,IB), FREQ_SP(IS), E_BT(IT_O,IS), E_OSF(IT_O,IS), DIFF1_ESP(IT_O,IS), DIFF2_ESP(IT_O,IS)
                  WRITE(49,'(4F20.8)') ZEIT_O(IT_O,IB), FREQ_SP(IS), E_BT(IT_O,IS), E_OSF(IT_O,IS)
                  !WRITE(50,'(4F20.8)') ZEIT_O(IT_O,IB), FREQ_SP(IS), E_OSF(IT_O,IS), E_BT(IT_O,IS)
                  !write(*,*) ib, ZEIT_O(IT_O,IB), FREQ_SP(IS), DIFF1_ESP(IT_O,IS), DIFF2_ESP(IT_O,IS)
                ELSE
                  DIFF1_ESP(IT_O,IS) =  0.0
                  DIFF2_ESP(IT_O,IS) =  0.0
                END IF
              END DO
              DO IS = 2, MSC
                IF((FREQ_SP(IS) .LE. FREQ_O(MSC_O(IB),IB))) THEN
                  ETMPS = ETMPS + E_BT(IT_O,IS)  * (FREQ_SP(IS)-FREQ_SP(IS-1))
                  ETMPO = ETMPO + E_OSF(IT_O,IS) * (FREQ_SP(IS)-FREQ_SP(IS-1)) 
                END IF
              END DO
              !WRITE(*,'(2I10,3F15.6)') IB, IT_O, ZEIT_O(IT_O,IB), 4*SQRT(ETMPS), 4*SQRT(ETMPO)
            ELSE
              DIFF1_ESP(IT_O,:) =  0.0
              DIFF2_ESP(IT_O,:) =  0.0
            END IF
          ELSE
            DIFF1_ESP(IT_O,:) =  0.0
            DIFF2_ESP(IT_O,:) =  0.0
          END IF ! Out of simulation period ...
          !IF (ZEIT_O(IT_O,IB) > ZEIT_SP(1) .AND. ZEIT_O(IT_O,IB) < ZEIT_SP(N_DT_SP)) THEN
          !  WRITE(48,'('//'F15.6,'//CHTMP//'F10.6'//')') ZEIT_O(IT_O,IB), (DIFF1_ESP(IT_O,IS), IS = 1, MSC)
          !END IF
        END DO
!
        !WRITE(*,*) IB, STATIONNAMES(IB)
        !WRITE(*,'(A6,5A16)') 'IT_O', 'FREQ', 'MEAN_BSP', 'MEAN_OSP', 'BIAS_SP', 'RMS_B'
        DO IS = 1, MSC
          IF((FREQ_SP(IS) .LE. FREQ_O(MSC_O(IB),IB))) THEN
            !MEAN_BSP(IB,IS) = SUM( E_BT(:,IS))/(N_DT_O(IB)-SUM(N_DT_ERR(IB,:)))
            !MEAN_OSP(IB,IS) = SUM(E_OSF(:,IS))/(N_DT_O(IB)-SUM(N_DT_ERR(IB,:)))
            !BIAS_SP(IB,IS)  = MEAN_BSP(IB,IS) - MEAN_OSP(IB,IS)
            !RMS_B(IB,IS)    = SQRT(SUM(DIFF2_ESP(:,IS)/((N_DT_O(IB)-SUM(N_DT_ERR(IB,:))))))
            MEAN_BSP(IB,IS) = SUM( E_BT(:,IS))/N_STAT(IB)!(N_DT_O(IB)-SUM(N_DT_ERR(IB,:)))
            MEAN_OSP(IB,IS) = SUM(E_OSF(:,IS))/N_STAT(IB)!(N_DT_O(IB)-SUM(N_DT_ERR(IB,:)))
            BIAS_SP(IB,IS)  = MEAN_BSP(IB,IS) - MEAN_OSP(IB,IS)
            RMS_B(IB,IS)    = SQRT(SUM(DIFF2_ESP(:,IS)/((N_STAT(IB)))))
            !WRITE(*,'(I10,5F15.4)') IS, FREQ_SP(IS), MEAN_BSP(IB,IS), MEAN_OSP(IB,IS), BIAS_SP(IB,IS), RMS_B(IB,IS)
          END IF
        END DO

        DO IT_O = 1, N_DT_O(IB)
          IF (ZEIT_O(IT_O,IB) > ZEIT_SP(1) .AND. ZEIT_O(IT_O,IB) < ZEIT_SP(N_DT_SP)) THEN
            IF (N_DT_ERR(IB,IT_O) .EQ. 0) THEN 
              DO IS = 1, MSC
                IF((FREQ_SP(IS) .LT. FREQ_O(MSC_O(IB),IB))) THEN
                  DIFF1_K(IT_O,IS) =  E_BT(IT_O,IS)  - MEAN_BSP(IB,IS)
                  DIFF2_K(IT_O,IS) = (E_BT(IT_O,IS)  - MEAN_BSP(IB,IS))**2.0
                  DIFF3_K(IT_O,IS) =  E_OSF(IT_O,IS) - MEAN_OSP(IB,IS)
                  DIFF4_K(IT_O,IS) = (E_OSF(IT_O,IS) - MEAN_OSP(IB,IS))**2.0
                ELSE
                  DIFF1_K(IT_O,IS) = 0.0
                  DIFF2_K(IT_O,IS) = 0.0
                  DIFF3_K(IT_O,IS) = 0.0
                  DIFF4_K(IT_O,IS) = 0.0
                END IF
              END DO
            ELSE
              DIFF1_K(IT_O,:) = 0.0
              DIFF2_K(IT_O,:) = 0.0
              DIFF3_K(IT_O,:) = 0.0
              DIFF4_K(IT_O,:) = 0.0
            END IF
          ELSE
            DIFF1_K(IT_O,:) = 0.0
            DIFF2_K(IT_O,:) = 0.0
            DIFF3_K(IT_O,:) = 0.0
            DIFF4_K(IT_O,:) = 0.0
          END IF
        END DO

        DO IS = 1, MSC
          IF((FREQ_SP(IS) .LE. FREQ_O(MSC_O(IB),IB))) THEN
            IF ((( SUM(DIFF2_K(:,IS))*SUM(DIFF4_K(:,IS)) )**0.5).GT.SMALL) THEN
              KORR_SP(IB,IS) = (SUM(DIFF1_K(:,IS)*DIFF3_K(:,IS))/((SUM(DIFF2_K(:,IS))*SUM(DIFF4_K(:,IS)))**0.5))**2.0
            ELSE
              KORR_SP(IB,IS) = 0.0
            END IF
          ELSE
            KORR_SP(IB,IS) = 0.0
!            WRITE(*,*) FREQ_SP(I), SUM( DIFF1_K(:,I) * DIFF3_K(:,I) ), (( SUM(DIFF2_K(:,I))  *  SUM(DIFF4_K(:,I)) )**0.5)            
          END IF        
        END DO
        
        DO IS = 1, MSC
          IF((FREQ_SP(IS) .LE. FREQ_O(MSC_O(IB),IB))) THEN
            WRITE(105,'(4F15.6)') FREQ_SP(IS), BIAS_SP(IB,IS), RMS_B(IB,IS), KORR_SP(IB,IS)
          END IF
        END DO
                
      END DO ! IB
           
      WRITE (101,*) 'KORRELATIONSKOEFFZIENTEN'
      WRITE (101,*)
      
      WRITE (101, '(32(A16))')  'FREQ', (BNAMES (IB), IB = 1, BUOYS)
      WRITE (101,*)
      
      DO IS = 1, MSC
        IF((FREQ_SP(IS) .LT. MINMAXFREQ_O)) THEN
          WRITE (101, '(F15.4, 31(F15.4))')  FREQ_SP(IS), (KORR_SP(IB,IS), IB = 1, BUOYS)
        END IF
      END DO
      
      WRITE (101,*) 'BIAS'
      WRITE (101,*)
      
      WRITE (101, '(32(A16))')  'FREQ', (BNAMES (IB), IB = 1, BUOYS)
      WRITE (101,*)
      
      DO IS = 1, MSC
        IF((FREQ_SP(IS) .LT. MINMAXFREQ_O)) THEN
          WRITE (101,'(F15.4, (40F15.4))')  FREQ_SP(IS), (BIAS_SP(IB,IS), IB = 1, BUOYS)
        END IF
      END DO
      
      WRITE (101,*) 'RMS'
      WRITE (101,*)
      
      WRITE (101, '(32(A16))')  'FREQ', (BNAMES (IB), IB = 1, BUOYS)
      WRITE (101,*)
      
      DO IS = 1, MSC
        IF((FREQ_SP(IS) .LT. MINMAXFREQ_O)) THEN
          WRITE (101,'(F15.4, (40F15.4))')  FREQ_SP(IS), (RMS_B(IB,IS), IB = 1, BUOYS)
        END IF
      END DO 
      
      WRITE (101,*)
      
      WRITE (101,*) 'FINAL STATISTICS FOR ALL BUOYS'
      WRITE (101,*)
            
      DO IS = 1, MSC
        IF(FREQ_SP(IS) .LT. MINMAXFREQ_O) THEN
          WRITE (101,'(4F15.6)') FREQ_SP(IS), SUM(BIAS_SP(:,IS))/IBUOYS, SUM(RMS_B(:,IS))/IBUOYS, SUM(KORR_SP(:,IS))/IBUOYS
        END IF
      END DO
      
      CLOSE (101) 
999   FORMAT(I0)
      END PROGRAM 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CT2MJD(STIME,XMJD)
         IMPLICIT NONE
         CHARACTER(LEN=15), INTENT(IN) :: STIME
         DOUBLE PRECISION, INTENT(INOUT) :: XMJD
         INTEGER :: IY, IM, ID, IH, IMIN, ISEC, IFLAG
!
! ... FORMAT IS YYYYMMDD.HHMMSS , LENGTH IS 15
!
         IFLAG = 1
         READ(STIME(1:4),*) IY
         READ(STIME(5:6),*) IM
         READ(STIME(7:8),*) ID
         READ(STIME(10:11),*) IH
         READ(STIME(12:13),*) IMIN
         READ(STIME(14:15),*) ISEC

         CALL MJDYMD(XMJD  , IY    , IM    , ID    , IH    ,   &
     &               IMIN  , ISEC  , IFLAG                    )

         RETURN
      END SUBROUTINE

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MJD2CT(XMJD,STIME)
         IMPLICIT NONE
         CHARACTER(LEN=15), INTENT(INOUT) :: STIME
         DOUBLE PRECISION, INTENT(IN) :: XMJD
         DOUBLE PRECISION :: TMJD
         INTEGER :: IY, IM, ID, IH, IMIN, ISEC, IFLAG

         IFLAG = 2
         TMJD = XMJD

         CALL MJDYMD(TMJD  , IY    , IM    , ID    , IH    ,   &
     &               IMIN  , ISEC  , IFLAG                      )

         WRITE(STIME(1:4),'(I4.4)') IY
         WRITE(STIME(5:6),'(I2.2)') IM
         WRITE(STIME(7:8),'(I2.2)') ID
         WRITE(STIME(9:9),'(A)') '.'
         WRITE(STIME(10:11),'(I2.2)') IH
         WRITE(STIME(12:13),'(I2.2)') IMIN
         WRITE(STIME(14:15),'(I2.2)') ISEC


         RETURN
      END SUBROUTINE

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CU2SEC(UNITT, DT)
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN) :: UNITT
         REAL*8, INTENT(INOUT) :: DT

         SELECT CASE (UNITT)
            CASE ('H', 'h', 'HR', 'hr')
               DT = DT * 3600.0
            CASE ('M', 'm', 'MIN', 'min')
               DT = DT * 60.0
            CASE ('S', 's', 'SEC', 'sec')
               DT = DT
            CASE DEFAULT
               WRITE(*,*) 'ERROR, UNIT = ', UNITT
               STOP
         END SELECT

         RETURN
      END SUBROUTINE

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MJDYMD(XMJD  , IY    , IM    , ID    , IH    ,  &
     &                  IMIN  , ISEC  , IFLAG                  )
!
! ... XMJD  : MODIFIED JULIAN DATE
! ... IY    : YEAR
! ... IM    : MONTH
! ... ID    : DAY
! ... IH    : HOUR
! ... IMIN  : MINUTE
! ... ISEC  : SECOND
! ... IFLAG : 1 -> YMDHMS TO MJD
!             2 -> MJD TO YMDHMS
! ... DATE MUST BE WITHIN THE YEARS MAR. 1, 1900 TO FEB. 28, 2100

         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
         PARAMETER ( XJD0 = 2400000.5D0 )
         PARAMETER ( HALF =       0.5D0 )
         INTEGER IMONTH(12)
         DATA IMONTH /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
!
! ... -----< YMDHMS TO MJD >-----
!
         IF (IFLAG.EQ.1) THEN

            Y = DFLOAT(IY - 1)

            IF (IM.GT.2) THEN
               M = IM
               Y = Y + 1
            ELSE
               M = IM + 12
            ENDIF

            XJD  = INT(365.25D0*Y) + INT(30.6001D0*(M+1)) - 15  &
     &           + 1720996.5D0     + ID
            XMJD = XJD - XJD0

            FSEC = DFLOAT(IH)*3600.D0 + DFLOAT(IMIN)*60.D0 + DFLOAT(ISEC)

            XMJD = XMJD + FSEC/86400.D0
!
! ... -----< MJD TO YMDHMS >-----
!
         ELSE IF (IFLAG.EQ.2) THEN

            MJD  = XMJD
            XJD  = DFLOAT(MJD) + XJD0
            C    = INT(XJD + HALF) + 1537
            ND   = INT((C - 122.1D0)/365.25D0 )
            E    = INT(365.25D0*ND)
            NF   = INT((C - E)/30.6001D0)

            IFR  = INT(XJD + HALF)
            FRC  = XJD + HALF - DFLOAT(IFR)
            ID   = C - E - INT(30.6001D0*NF) + FRC
            IM   = NF - 1 - 12*INT(NF/14)
            IY   = ND - 4715 - INT((7+IM)/10)

            SEC  = (XMJD-DFLOAT(MJD))*86400.D0
            ISEC = SEC
            IF ((SEC-ISEC).GT.0.5D0) ISEC = ISEC + 1
            IH   = ISEC/3600
            IMIN = (ISEC - IH*3600)/60
            ISEC = ISEC - IH*3600 - IMIN*60
!
! ... set 24:00 to 00:00
!
            IF (MOD(IY,4) == 0) IMONTH(2) = 29

            IF (IH == 24) THEN
               IH = 0
               ID = ID + 1
               IF (ID > IMONTH(IM)) THEN
                  ID = ID - IMONTH(IM)
                  IM = IM + 1
                  IF (IM > 12) THEN
                     IM = IM - 12
                     IY = IY + 1
                  END IF
               END IF
            END IF
         ELSE

           PRINT*,'!!! ERROR IN <MJDYMD>. IFLAG SHOULD BE 1 OR 2.'
           STOP

         ENDIF

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************      
      SUBROUTINE INTER_S (Y1,Y2,DX,DIFF_DX,YINTER)
      IMPLICIT NONE

      REAL*8, INTENT(OUT) :: YINTER
      REAL*8, INTENT(IN)  :: Y1, Y2, DX, DIFF_DX

      YINTER = Y1 + DIFF_DX * (Y2-Y1) / DX
      IF (DIFF_DX == 0.0) YINTER = Y1

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_PARAMETER_LOC(ACLOC,ACLOC2D,CURTXYLOC,DEPLOC,WKLOC,ISMAX,HS,TM01,TM02,KLM,WLM)
         USE STAT_POOL
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: ISMAX
         REAL, INTENT(IN)    :: ACLOC(MSC,MDC),ACLOC2D(MSC,MDC)
         REAL, INTENT(IN)    :: WKLOC(MSC), DEPLOC
         REAL, INTENT(IN)    :: CURTXYLOC(2)


         REAL*8, INTENT(OUT)   :: HS,TM01,TM02,KLM,WLM

         INTEGER               :: ID, IS

         REAL*8                :: Y(MSC), CURMAG
         REAL*8                :: OMEG2,OMEG,EAD,UXD,UYD,ETOT,EFTOT,APTOT,EPTOT
         REAL*8                :: SKK, CKTAIL, ETOT1, SIG22, EKTOT, CETAIL
         REAL*8                :: dintspec, dintspec_y, tmp(msc),actmp(msc)
!
! total energy ...
!
         ETOT = 0.
         do id = 1, mdc
           tmp = acloc(:,id) * spsig
           do is = 2, ismax
             ETOT = ETOT + 0.5*(tmp(is)+tmp(is-1))*(spsig(is)-spsig(is-1))*ddir
           end do
         end do

         HS = 4*SQRT(ETOT)

         APTOT = 0.
         EPTOT = 0.

         DO ID = 1, MDC
           DO IS = 1, ISMAX
             APTOT = APTOT + SPSIG(IS)     * ACLOC2D(IS,ID)
             EPTOT = EPTOT + SPSIG(IS)**2. * ACLOC2D(IS,ID)
           ENDDO
         ENDDO

         APTOT = APTOT * FRINTF
         EPTOT = EPTOT * FRINTF

         IF (EPTOT .GT. SMALL) THEN
            TM01 = 2.*PI * APTOT / EPTOT
         ELSE
            TM01 = 0.
         END IF

         ETOT = 0.
         EFTOT = 0.
         DO ID=1, MDC
            DO IS = 1, ISMAX
              EAD   = SPSIG(IS) * ACLOC2D(IS,ID) * FRINTF
              OMEG2 = SPSIG(IS)**2
              ETOT  = ETOT + EAD
              EFTOT = EFTOT + EAD * OMEG2
            ENDDO
         ENDDO

         IF (EFTOT .GT. SMALL .AND. ETOT .GT. SMALL) THEN
           TM02 = PI2 * SQRT(ETOT/EFTOT)
         ELSE
           TM02 = 0.
         END IF

         ETOT1 = 0.
         EKTOT = 0.
!
         DO IS = 1, ISMAX
           SIG22 = SPSIG(IS)**2.
           SKK  = SIG22 * WKLOC(IS)
           DO ID = 1, MDC
             ETOT1 = ETOT1 + SIG22 * ACLOC2D(IS,ID)
             EKTOT = EKTOT + SKK * ACLOC2D(IS,ID)
           ENDDO
         ENDDO

         ETOT1 = FRINTF * ETOT1
         EKTOT = FRINTF * EKTOT

         IF (ETOT1.GT.SMALL.AND.EKTOT.GT.SMALL) THEN
            WLM = PI2 * (ETOT1/EKTOT)
            KLM  = PI2/WLM
         ELSE
            KLM  = 10.
            WLM = 0.
         ENDIF

       RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************    
      SUBROUTINE MEAN_PARAMETER_LOC_ACTION(ACLOC,CURTXYLOC,DEPLOC,WKLOC,ISMAX,HS,TM01,TM02,KLM,WLM)
         USE STAT_POOL
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: ISMAX
         REAL, INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL, INTENT(IN)    :: WKLOC(MSC), DEPLOC
         REAL, INTENT(IN)    :: CURTXYLOC(2)


         REAL*8, INTENT(OUT)   :: HS,TM01,TM02,KLM,WLM

         INTEGER               :: ID, IS

         REAL*8                :: Y(MSC)
         REAL*8                :: OMEG2,OMEG,EAD,UXD,UYD,ETOT,EFTOT,APTOT,EPTOT
         REAL*8                :: SKK, CKTAIL, ETOT1, SIG22, EKTOT, CETAIL
         REAL*8                :: dintspec, dintspec_y, tmp(msc),actmp(msc)
!
! total energy ...
!
         ETOT = 0.
         do id = 1, mdc
           tmp = acloc(:,id) * spsig
           do is = 2, ismax
             ETOT = ETOT + 0.5*(tmp(is)+tmp(is-1))*(spsig(is)-spsig(is-1))*ddir
           end do
         end do

         HS = 4*SQRT(ETOT)

         APTOT = 0.
         EPTOT = 0.

         DO ID = 1, MDC
           DO IS = 1, ISMAX
             APTOT = APTOT + SPSIG(IS) * ACLOC(IS,ID)
             EPTOT = EPTOT + SPSIG(IS)**2. * ACLOC(IS,ID)
           ENDDO
         ENDDO

         APTOT = APTOT * FRINTF
         EPTOT = EPTOT * FRINTF

         IF (EPTOT .GT. SMALL) THEN
            TM01 = 2.*PI * APTOT / EPTOT
         ELSE
            TM01 = 0.
         END IF


         ETOT  = 0.
         EFTOT = 0.
         DO ID=1, MDC
            IF (CURTXYLOC(1)**2.+CURTXYLOC(2)**2.GT.SMALL) THEN
              UXD  = CURTXYLOC(1)*COSTH(ID) + CURTXYLOC(2)*SINTH(ID)
            ENDIF
            DO IS = 1, ISMAX
              EAD  = SPSIG(IS)**2 * ACLOC(IS,ID) * FRINTF
              IF (CURTXYLOC(1)**2.+CURTXYLOC(2)**2.GT.SMALL) THEN
                OMEG  = SPSIG(IS) + WKLOC(IS) * UXD
                OMEG2 = OMEG**2
              ELSE
               OMEG2 = SPSIG(IS)**2
              ENDIF
              ETOT  = ETOT + EAD
              EFTOT = EFTOT + EAD * OMEG2
            ENDDO
         ENDDO
         IF (EFTOT .GT. SMALL) THEN
           TM02 = PI2 * SQRT(ETOT/EFTOT)
         ELSE
           TM02 = 0.
         END IF

         ETOT1 = 0.
         EKTOT = 0.
!
         DO IS = 1, ISMAX
           SIG22 = SPSIG(IS)**2.
           SKK  = SIG22 * WKLOC(IS)
           DO ID = 1, MDC
             ETOT1 = ETOT1 + SIG22 * ACLOC(IS,ID)
             EKTOT = EKTOT + SKK * ACLOC(IS,ID)
           ENDDO
         ENDDO

         ETOT1 = FRINTF * ETOT1
         EKTOT = FRINTF * EKTOT

         IF (ETOT1.GT.SMALL.AND.EKTOT.GT.SMALL) THEN
            WLM = PI2 * (ETOT1/EKTOT)
            KLM  = PI2/WLM
         ELSE
            KLM  = 10.
            WLM = 0.
         ENDIF
       RETURN
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
        Print *, TRIM(ErrMsg)
        STOP
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************    

