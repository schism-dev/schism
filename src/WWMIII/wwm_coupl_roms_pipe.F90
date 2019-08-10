#include "wwm_functions.h"
#if !defined ROMS_WWM_PGMCL_COUPLING && !defined MODEL_COUPLING_ATM_WAV && !defined MODEL_COUPLING_OCN_WAV
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_PIPES_ROMS()
      USE DATAPOOL
      IMPLICIT NONE
!
! open pipe data files for coupling
!
      LSEWL = .TRUE.
      LSECU = .TRUE.
      WRITE(DBG%FHNDL,'("+TRACE...",A)') 'OPEN PIPE ROMS'
      FLUSH(DBG%FHNDL)
!     Pipes that are read by the wave model
      OPEN(1000,file='pipe/ExchRW'  ,form='unformatted', action='read')
      WRITE(DBG%FHNDL,*) 'WWM: open pipe ExchImport'
      FLUSH(DBG%FHNDL)
!     Pipes that are written by the wave modell
      OPEN(101 ,file='pipe/ExchWR' ,form='unformatted', action='write')
      WRITE(DBG%FHNDL,*) 'WWM: open pipe ExchExport'
      FLUSH(DBG%FHNDL)
      WRITE(DBG%FHNDL,'("+TRACE...",A)') 'END OPEN PIPE ROMS'
      FLUSH(DBG%FHNDL)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE TERMINATE_PIPES_ROMS()
      USE DATAPOOL
      IMPLICIT NONE
      close(1000)
      close(101)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PIPE_ROMS_IN(K,IFILE,IT)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: K,IFILE,IT
      INTEGER              :: IP
# ifdef WWM_MPI
      REAL(rkind), allocatable :: WINDXY_TOT(:,:), CURTXY_TOT(:,:), WATLEV_TOT(:)
      real(rkind), allocatable :: rbuf_real(:)
      integer idx, iProc
# endif
        LCALC=.TRUE.
      IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN
        WATLEVOLD=WATLEV
        LCALC=.TRUE.
        WRITE(DBG%FHNDL,'("+TRACE...",A)') 'READING PIPE'
        FLUSH(DBG%FHNDL)
# ifndef WWM_MPI
        DO IP = 1, MNP
          READ(1000) WINDXY(IP,1), WINDXY(IP,2), CURTXY(IP,1), CURTXY(IP,2), WATLEV(IP)
        END DO
# else
        allocate(WINDXY_TOT(np_global,2), CURTXY_TOT(np_global,2), WATLEV_TOT(np_global), rbuf_real(np_global*5), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate err')
        IF (myrank.eq.0) THEN
          DO IP = 1, np_global
            READ(1000) WINDXY_TOT(IP,1), WINDXY_TOT(IP,2), CURTXY_TOT(IP,1), CURTXY_TOT(IP,2), WATLEV_TOT(IP)
          END DO
          DO IP=1,np_global
            idx=idx+1
            rbuf_real(idx)=WINDXY_TOT(IP,1)
            idx=idx+1
            rbuf_real(idx)=WINDXY_TOT(IP,2)
            idx=idx+1
            rbuf_real(idx)=CURTXY_TOT(IP,1)
            idx=idx+1
            rbuf_real(idx)=CURTXY_TOT(IP,2)
            idx=idx+1
            rbuf_real(idx)=WATLEV_TOT(IP)
          END DO
          DO iProc=2,nproc
            CALL MPI_SEND(rbuf_real,np_global*5,MPI_REAL8, iProc-1, 196, comm, ierr)
          END DO
        ELSE
          CALL MPI_RECV(rbuf_real,np_global*5,MPI_REAL8, 0, 196, comm, istatus, ierr)
          idx=0
          DO IP=1,np_global
            idx=idx+1
            WINDXY_TOT(IP,1)=rbuf_real(idx)
            idx=idx+1
            WINDXY_TOT(IP,2)=rbuf_real(idx)
            idx=idx+1
            CURTXY_TOT(IP,1)=rbuf_real(idx)
            idx=idx+1
            CURTXY_TOT(IP,2)=rbuf_real(idx)
            idx=idx+1
            WATLEV_TOT(IP)=rbuf_real(idx)
          END DO
        END IF
        DO IP = 1, MNP
          WINDXY(IP,:)=WINDXY_TOT(iplg(IP),:)
          CURTXY(IP,:)=CURTXY_TOT(iplg(IP),:)
          WATLEV(IP)=WATLEV_TOT(iplg(IP))
        END DO
        deallocate(rbuf_real, WINDXY_TOT, CURTXY_TOT, WATLEV_TOT)
# endif
        DEPDT = (WATLEV - WATLEVOLD) / MAIN%DTCOUP
        WRITE(DBG%FHNDL,'("+TRACE...",A)') 'END READING PIPE'
        FLUSH(DBG%FHNDL)
      END IF
      IF (K == 1) CALL INITIAL_CONDITION
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PIPE_ROMS_OUT(K)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: K
      INTEGER              :: IP
      REAL(rkind)          :: ACLOC(MSC,MDC)
      REAL(rkind)          :: HS,WLM,LPP,FPP,CPP,BOTEXPER
      REAL(rkind)          :: UBOT,TM01,TM10
      REAL(rkind)          :: TMBOT, KPP,DM,DSPR,ORBITAL,ETOTS,ETOTC,WNPP,TPP,CGPP
      REAL(rkind)          :: PEAKDSPR, PEAKDM, HSWE, HSLIM, TM02, KLM, DPEAK
      REAL(rkind)          :: TPPD,KPPD,CGPD,CPPD
# ifdef WWM_MPI
      REAL(rkind), allocatable :: OUTT(:,:), OUTT_TOT(:,:)
      REAL(rkind)    :: TP
# endif
      IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN
# ifndef WWM_MPI
        DO IP = 1, MNP
          ACLOC = AC2(:,:,IP)
          CALL MEAN_PARAMETER(IP,ACLOC,MSC,HS,TM01,TM02,TM10,KLM,WLM)
          CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT,'PIPE_ROMS_OUT 1')
          CALL MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,MSC,ETOTS,ETOTC,DM,DSPR)
          CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
! HS, HSWE, HSLIM  ! - Significant wave height (m) -- HS
! DM  ! - Wave direction (degrees)
! TPP ! - Surface wave relative peak period (s) -- TP
! WLM, KME, SME ! - Average Wave Length [m] - LME
! ORBITAL(IP) ! - Wave bottom orbital velocity (m/s)
! TMBOT ! - Bottom wave period (s)
! DISSIPATION(IP) ! - Wave energy dissipation (W/m2)
! QBLOCAL(IP) ! - Percent of breakig waves (nondimensional)
! DSPR ! - directional spreading
! PEAKDSPR ! - peak directional spreading
! PEAKDM ! - Peak direction
!
!AR: what for you need this HSWE and HSLIM and so on ... what is exactly SME in your definition ...
!AR: I have deleted them ...
!
          HSWE=0
          HSLIM=0
          WRITE(101)  HS, HSWE,                                        &
     &                   HSLIM, DM,                                    &
     &                   TPP, WLM,                                     &
     &                   KLM, TM01,                                    &
     &                   ORBITAL, TMBOT,                               &
     &                   DISSIPATION(IP), QBLOCAL(IP),                 &
     &                   DSPR, PEAKDSPR,                               &
     &                   PEAKDM, TM02
          FLUSH(101)
        END DO
# else
        allocate(OUTT(np_global,16), OUTT_TOT(np_global,16), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate err')
        OUTT=0
        DO IP = 1, MNP
          ACLOC = AC2(:,:,IP)
          CALL MEAN_PARAMETER(IP,ACLOC,MSC,HS,TM01,TM02,TM10,KLM,WLM)
          CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT,'PIPE_ROMS_OUT 2')
          CALL MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,MSC,ETOTS,ETOTC,DM,DSPR)
          CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
          HSWE=0
          HSLIM=0
          OUTT(iplg(IP), 1)=HS
          OUTT(iplg(IP), 2)=HSWE
          OUTT(iplg(IP), 3)=HSLIM
          OUTT(iplg(IP), 4)=DM
          OUTT(iplg(IP), 5)=TPP
          OUTT(iplg(IP), 6)=WLM
          OUTT(iplg(IP), 7)=KLM
          OUTT(iplg(IP), 8)=TM01
          OUTT(iplg(IP), 9)=ORBITAL
          OUTT(iplg(IP),10)=TMBOT
          OUTT(iplg(IP),11)=DISSIPATION(IP)
          OUTT(iplg(IP),12)=QBLOCAL(IP)
          OUTT(iplg(IP),13)=DSPR
          OUTT(iplg(IP),14)=PEAKDSPR
          OUTT(iplg(IP),15)=PEAKDM
          OUTT(iplg(IP),16)=TM02
        END DO
        call mpi_reduce(OUTT,OUTT_TOT,NP_GLOBAL*16,rtype,MPI_SUM,0,comm,ierr)
        IF (myrank.eq.0) THEN
          DO IP=1,NP_GLOBAL
            OUTT_TOT(IP,:)=OUTT_TOT(IP,:)/nwild_gb(IP)
            WRITE(101) OUTT_TOT(IP, 1), OUTT_TOT(IP, 2),                 &
     &                    OUTT_TOT(IP, 3), OUTT_TOT(IP, 4),              &
     &                    OUTT_TOT(IP, 5), OUTT_TOT(IP, 6),              &
     &                    OUTT_TOT(IP, 7), OUTT_TOT(IP, 8),              &
     &                    OUTT_TOT(IP, 9), OUTT_TOT(IP,10),              &
     &                    OUTT_TOT(IP,11), OUTT_TOT(IP,12),              &
     &                    OUTT_TOT(IP,13), OUTT_TOT(IP,14),              &
     &                    OUTT_TOT(IP,15), OUTT_TOT(IP,16)
          END DO
        END IF
        deallocate(OUTT, OUTT_TOT)
# endif
      END IF
      WRITE(DBG%FHNDL,*) 'export WWM: ending of writing data'
      FLUSH(DBG%FHNDL)
      END SUBROUTINE
#endif
