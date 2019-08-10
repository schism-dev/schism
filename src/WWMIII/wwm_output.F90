#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GENERAL_OUTPUT
      USE WWM_HOTFILE_MOD
      USE DATAPOOL
      IMPLICIT NONE
      WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4,L5)') 'WRITING OUTPUT INTERNAL TIME', RTIME, MAIN%TMJD, OUT_HISTORY%TMJD-1.E-8, OUT_HISTORY%EMJD, (MAIN%TMJD .GE. OUT_HISTORY%TMJD-1.E-8) .AND. (MAIN%TMJD .LE. OUT_HISTORY%EMJD)
      !
      ! The history output
      !
      IF ( (MAIN%TMJD .GE. OUT_HISTORY%TMJD-1.E-8) .AND. (MAIN%TMJD .LE. OUT_HISTORY%EMJD)) THEN
        WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'WRITING OUTPUT INTERNAL TIME', RTIME, MAIN%TMJD, OUT_HISTORY%TMJD-1.E-8, OUT_HISTORY%EMJD
        CALL OUTPUT_HISTORY(RTIME*DAY2SEC,.FALSE.)
        OUT_HISTORY%TMJD = OUT_HISTORY%TMJD + OUT_HISTORY%DELT*SEC2DAY
      END IF
      !
      ! The station output
      !
      IF ( (MAIN%TMJD .GE. OUT_STATION%TMJD-1.E-8) .AND. (MAIN%TMJD .LE. OUT_STATION%EMJD)) THEN
        WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)')  'WRITING OUTPUT INTERNAL TIME', RTIME, MAIN%TMJD, OUT_STATION%TMJD-1.E-8, OUT_STATION%EMJD
        CALL OUTPUT_STATION(RTIME*DAY2SEC,.FALSE.)
        OUT_STATION%TMJD = OUT_STATION%TMJD + OUT_STATION%DELT*SEC2DAY
      END IF
      !
      ! The hotfile output
      !
      IF (LHOTF) THEN
        IF ( (MAIN%TMJD .GE. HOTF%TMJD-1.E-8) .AND. (MAIN%TMJD .LE. HOTF%EMJD)) THEN
          WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'WRITING HOTFILE INTERNAL TIME', RTIME
          FLUSH(STAT%FHNDL)
          CALL OUTPUT_HOTFILE
          HOTF%TMJD = HOTF%TMJD + HOTF%DELT*SEC2DAY
        END IF
      END IF
      !
      ! The wavewatch III exports
      !
      WRITE(STAT%FHNDL,*) 'Before LEXPORT_BOUC_WW3'
      FLUSH(STAT%FHNDL)
      IF (LEXPORT_BOUC_WW3) THEN
        WRITE(STAT%FHNDL,*) 'Before time test'
        WRITE(STAT%FHNDL,*) 'MAIN%TMJD=', MAIN%TMJD
        WRITE(STAT%FHNDL,*) 'OUT_BOUC_WW3%TMJD=', OUT_BOUC_WW3%TMJD
        WRITE(STAT%FHNDL,*) 'OUT_BOUC_WW3%EMJD=', OUT_BOUC_WW3%EMJD
        FLUSH(STAT%FHNDL)
        IF ( (MAIN%TMJD .GE. OUT_BOUC_WW3%TMJD-1.E-8) .AND. (MAIN%TMJD .LE. OUT_BOUC_WW3%EMJD)) THEN
          WRITE(STAT%FHNDL,*) 'After time test'
          FLUSH(STAT%FHNDL)
!          CALL EXPORT_BOUC_WW3_FORMAT
          OUT_BOUC_WW3%TMJD = OUT_BOUC_WW3%TMJD + OUT_BOUC_WW3%DELT*SEC2DAY
        END IF
      END IF
      WRITE(STAT%FHNDL,*) 'Before LEXPORT_WIND_WW3'
      FLUSH(STAT%FHNDL)
      IF (LEXPORT_WIND_WW3) THEN
        WRITE(STAT%FHNDL,*) 'Before time test'
        FLUSH(STAT%FHNDL)
        IF ( (MAIN%TMJD .GE. OUT_WIND_WW3%TMJD-1.E-8) .AND. (MAIN%TMJD .LE. OUT_WIND_WW3%EMJD)) THEN
          WRITE(STAT%FHNDL,*) 'After time test'
          FLUSH(STAT%FHNDL)
!          CALL EXPORT_WIND_WW3_FORMAT
          OUT_WIND_WW3%TMJD = OUT_WIND_WW3%TMJD + OUT_WIND_WW3%DELT*SEC2DAY
        END IF
      END IF
      WRITE(STAT%FHNDL,*) 'Before LEXPORT_CURR_WW3'
      FLUSH(STAT%FHNDL)
      IF (LEXPORT_CURR_WW3) THEN
        WRITE(STAT%FHNDL,*) 'Before time test'
        FLUSH(STAT%FHNDL)
        IF ( (MAIN%TMJD .GE. OUT_CURR_WW3%TMJD-1.E-8) .AND. (MAIN%TMJD .LE. OUT_CURR_WW3%EMJD)) THEN
          WRITE(STAT%FHNDL,*) 'After time test'
          FLUSH(STAT%FHNDL)
!          CALL EXPORT_CURR_WW3_FORMAT
          OUT_CURR_WW3%TMJD = OUT_CURR_WW3%TMJD + OUT_CURR_WW3%DELT*SEC2DAY
        END IF
      END IF
      WRITE(STAT%FHNDL,*) 'Before LEXPORT_WALV_WW3'
      FLUSH(STAT%FHNDL)
      IF (LEXPORT_WALV_WW3) THEN
        WRITE(STAT%FHNDL,*) 'Before time test'
        FLUSH(STAT%FHNDL)
        IF ( (MAIN%TMJD .GE. OUT_WALV_WW3%TMJD-1.E-8) .AND. (MAIN%TMJD .LE. OUT_WALV_WW3%EMJD)) THEN
          WRITE(STAT%FHNDL,*) 'After time test'
          FLUSH(STAT%FHNDL)
!          CALL EXPORT_WALV_WW3_FORMAT
          OUT_WALV_WW3%TMJD = OUT_WALV_WW3%TMJD + OUT_WALV_WW3%DELT*SEC2DAY
        END IF
      END IF
      !
      ! The boundary output
      !
      IF (BOUC_NETCDF_OUT_PARAM .or. BOUC_NETCDF_OUT_SPECTRA) THEN
#ifdef NCDF
        IF ( (MAIN%TMJD .GE. OUT_BOUC%TMJD-1.E-8) .AND. (MAIN%TMJD .LE. OUT_BOUC%EMJD)) THEN
          WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'WRITING OUTPUT INTERNAL TIME', RTIME, MAIN%TMJD, OUT_BOUC%TMJD-1.E-8, OUT_BOUC%EMJD
          CALL WRITE_NETCDF_BOUNDARY
          OUT_BOUC%TMJD = OUT_BOUC%TMJD + OUT_BOUC%DELT*SEC2DAY
        END IF
#else
        CALL WWM_ABORT('Need netcdf for the boundary output')
#endif
      END IF
      !
      ! The nesting
      !
      IF (L_NESTING) THEN
#ifdef NCDF
        CALL DO_NESTING_OPERATION
#else
        CALL WWM_ABORT('Need netcdf for the nesting output')
#endif
      END IF
      !
      WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH GENERAL_OUTPUT'
      FLUSH(STAT%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WWM_OUTPUT( TIME, LINIT_OUTPUT )
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(IN)    :: TIME
      LOGICAL, INTENT(IN) :: LINIT_OUTPUT
      CALL OUTPUT_HISTORY( TIME, LINIT_OUTPUT )
      CALL OUTPUT_STATION( TIME, LINIT_OUTPUT )
      OUT_HISTORY%TMJD = OUT_HISTORY%TMJD + OUT_HISTORY%DELT*SEC2DAY
      OUT_STATION%TMJD = OUT_STATION%TMJD + OUT_STATION%DELT*SEC2DAY

      WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH WWM OUTPUT'
      FLUSH(STAT%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OUTPUT_HISTORY( TIME, LINIT_OUTPUT )
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(IN)    :: TIME
      LOGICAL, INTENT(IN) :: LINIT_OUTPUT
      SELECT CASE (VAROUT_HISTORY%IOUTP)
        CASE (0)
          ! Do nothing ...
        CASE (1)
          CALL OUTPUT_HISTORY_XFN( TIME, LINIT_OUTPUT )
        CASE (2)
#ifdef NCDF
          CALL OUTPUT_HISTORY_NC
#else
          CALL WWM_ABORT('For History in netcdf, need netcdf!')
#endif
#ifdef DARKO
        CASE (3)
          CALL OUTPUT_HISTORY_SHP( TIME )
#endif
        CASE DEFAULT
          WRITE(DBG%FHNDL,*) 'IOUTP=', VAROUT_HISTORY%IOUTP
          CALL WWM_ABORT('WRONG NO OUTPUT SPECIFIED')
      END SELECT
      WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH OUTPUT_HISTORY'
      FLUSH(STAT%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OUTPUT_STATION( TIME, LINIT_OUTPUT )
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(IN)    :: TIME
      LOGICAL, INTENT(IN) :: LINIT_OUTPUT
      CHARACTER(LEN=15)   :: CTIME
      CALL MJD2CT(MAIN%TMJD, CTIME)
      IF ((DIMMODE .GT. 1) .and. LOUTS) THEN
        WRITE(STAT%FHNDL,*) 'WRITING STATION OUTPUT'
        SELECT CASE (VAROUT_STATION%IOUTP)
          CASE (0)
            ! Do nothing ...
          CASE (1)
            CALL OUTPUT_STE(CTIME, LINIT_OUTPUT)
          CASE (2)
#ifdef NCDF
            CALL OUTPUT_STATION_NC
#else
            CALL WWM_ABORT('STATION_NC: Need to compile with netcdf')
#endif
          CASE DEFAULT
            CALL WWM_ABORT('WRONG NO STATION OUTPUT SPECIFIED')
        END SELECT
      END IF
      WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH OUTPUT_HISTORY STATION'        
      FLUSH(STAT%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OUTPUT_HISTORY_XFN( TIME, LINIT_OUTPUT )
!
!     XFN TYPE OUTPUT
!
         USE DATAPOOL
         IMPLICIT NONE
         REAL(rkind), INTENT(IN)   :: TIME
         ! Yes we really want kind=4 variables here. The xfn tools can read kind=4 only
         LOGICAL, INTENT(IN)       :: LINIT_OUTPUT

         INTEGER                   :: IP
         LOGICAL                   :: DoAirSea
#ifdef MPI_PARALL_GRID
         REAL(rkind)               :: OUTT_GLOBAL(NP_GLOBAL,OUTVARS)
         REAL(rkind)               :: OUTT(NP_GLOBAL,OUTVARS)
         REAL(rkind)               :: CURR_GLOBAL(NP_GLOBAL,CURRVARS)
         REAL(rkind)               :: CURR(NP_GLOBAL,CURRVARS)
         REAL(rkind)               :: FORCE_GLOBAL(NP_GLOBAL,2)
         REAL(rkind)               :: FORCE(NP_GLOBAL,2)
         REAL(rkind)               :: WIND_GLOBAL(NP_GLOBAL,WINDVARS)
         REAL(rkind)               :: WIND(NP_GLOBAL,WINDVARS)
         REAL(rkind)               :: ITER_GLOBAL(NP_GLOBAL), ITER_LOCAL(MNP)
         REAL(rkind)               :: ITERT(NP_GLOBAL)
#else
         REAL(rkind)               :: OUTT(MNP,OUTVARS)
         REAL(rkind)               :: CURR(MNP,CURRVARS)
         REAL(rkind)               :: WIND(MNP,WINDVARS)
#endif
         REAL(rkind)               :: ACLOC(MSC,MDC)
         REAL(rkind)               :: OUTPARS(OUTVARS)
         REAL(rkind)               :: CURRPARS(CURRVARS)
         REAL(rkind)               :: WINDPARS(WINDVARS)

         CHARACTER(LEN=15)  :: CTIME

         CALL MJD2CT(MAIN%TMJD, CTIME)
         DoAirSea=.FALSE.
#ifdef MPI_PARALL_GRID
         OUTT_GLOBAL  = zero 
         OUTT         = zero
         CURR_GLOBAL  = zero
         CURR         = zero
         WIND_GLOBAL  = zero
         WIND         = zero
         FORCE_GLOBAL = zero
         FORCE        = zero
#else
         OUTT = zero 
         CURR = zero 
         WIND = zero 
#endif
         ACLOC    = zero 
         OUTPARS  = zero 
         CURRPARS = zero 
         WINDPARS = zero 

#ifdef MPI_PARALL_GRID
         IF (LQSTEA .AND. LCHKCONV) ITER_LOCAL = MyREAL(IP_IS_STEADY)
#endif

#ifdef MPI_PARALL_GRID
         DO IP = 1, MNP
            IF (DEP(IP) .GT. DMIN) THEN
              ACLOC(:,:) = AC2(:,:,IP)
              CALL INTPAR(IP, MSC, ACLOC, OUTPARS)
              CALL CURRPAR(IP, CURRPARS)
              CALL WINDPAR(IP, WINDPARS)
              IF (LMONO_OUT) OUTPARS(1) = OUTPARS(1) / SQRT(2.)
            ELSE
              OUTPARS     = zero 
              CURRPARS    = zero 
              WINDPARS    = zero 
            END IF
            OUTT(iplg(IP),:)   = OUTPARS(:)
            CURR(iplg(IP),:)   = CURRPARS(:)
            WIND(iplg(IP),:)   = WINDPARS(:)
            FORCE(iplg(IP),:)  = FORCEXY(IP,:)
            IF (LQSTEA) ITERT(iplg(IP))  = ITER_LOCAL(IP)
         END DO

         call mpi_reduce(OUTT,OUTT_GLOBAL,NP_GLOBAL*OUTVARS,rtype,MPI_SUM,0,comm,ierr)
         call mpi_reduce(CURR,CURR_GLOBAL,NP_GLOBAL*CURRVARS,rtype,MPI_SUM,0,comm,ierr)
         call mpi_reduce(WIND,WIND_GLOBAL,NP_GLOBAL*WINDVARS,rtype,MPI_SUM,0,comm,ierr)
         call mpi_reduce(FORCE,FORCE_GLOBAL,NP_GLOBAL*2,rtype,MPI_SUM,0,comm,ierr)
         IF (LQSTEA  .AND. LCHKCONV) call mpi_reduce(ITERT,ITER_GLOBAL,NP_GLOBAL,rtype,MPI_SUM,0,comm,ierr)

         if(myrank==0) then
           do IP=1,NP_GLOBAL
             OUTT_GLOBAL(IP,:) = OUTT_GLOBAL(IP,:)*nwild_gb(IP)
             CURR_GLOBAL(IP,:) = CURR_GLOBAL(IP,:)*nwild_gb(IP)
             WIND_GLOBAL(IP,:) = WIND_GLOBAL(IP,:)*nwild_gb(IP)
             FORCE_GLOBAL(IP,:) = FORCE_GLOBAL(IP,:)*nwild_gb(IP)
             IF (LQSTEA .AND. LCHKCONV) ITER_GLOBAL(IP)  =ITER_GLOBAL(IP)  *nwild_gb(IP)
           enddo !IP
         endif !myrank

         IF (myrank == 0) THEN
           IF (LINIT_OUTPUT) THEN
             OPEN(OUT%FHNDL+1, FILE  = 'ergzusw.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+2, FILE  = 'erguvh.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+3, FILE  = 'ergwind.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+4, FILE  = 'ergtm02.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+5, FILE  = 'erguvd.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+6, FILE  = 'ergufric.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+7, FILE  = 'ergtau.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+8, FILE  = 'ergiter.bin'  , FORM = 'UNFORMATTED')
             IF (DoAirSea) THEN
               OPEN(OUT%FHNDL+9, FILE  = 'airsea.dat'  , FORM = 'FORMATTED')
             END IF
           END IF
           WRITE(OUT%FHNDL+1)  SNGL(TIME)
           !WRITE(OUT%FHNDL+1)  (SNGL(FORCE_GLOBAL(IP,1)), SNGL(FORCE_GLOBAL(IP,2)), SNGL(OUTT_GLOBAL(IP,1))  , IP = 1, NP_GLOBAL)
           WRITE(OUT%FHNDL+1)  (SNGL(OUTT_GLOBAL(IP,7)), SNGL(OUTT_GLOBAL(IP,8)), SNGL(OUTT_GLOBAL(IP,1))  , IP = 1, NP_GLOBAL)
           CALL FLUSH(OUT%FHNDL+1)
           WRITE(OUT%FHNDL+2)  SNGL(TIME)
           !WRITE(OUT%FHNDL+2)  (SNGL(CURR_GLOBAL(IP,1)), SNGL(CURR_GLOBAL(IP,2)), SNGL(ZETA_SETUP(IP)), IP = 1, NP_GLOBAL)
           WRITE(OUT%FHNDL+2)  (SNGL(CURR_GLOBAL(IP,1)), SNGL(CURR_GLOBAL(IP,2)), SNGL(CURR_GLOBAL(IP,3)), IP = 1, NP_GLOBAL)
           FLUSH(OUT%FHNDL+2)
           WRITE(OUT%FHNDL+3)  SNGL(TIME)
           WRITE(OUT%FHNDL+3)  (SNGL(WIND_GLOBAL(IP,1)), SNGL(WIND_GLOBAL(IP,2)), SNGL(OUTT_GLOBAL(IP,10))  , IP = 1, NP_GLOBAL)
           FLUSH(OUT%FHNDL+3)
           WRITE(OUT%FHNDL+4)  SNGL(TIME)
           WRITE(OUT%FHNDL+4)  (SNGL(OUTT_GLOBAL(IP,1)), SNGL(OUTT_GLOBAL(IP,2)), SNGL(OUTT_GLOBAL(IP,3))  , IP = 1, NP_GLOBAL)
           FLUSH(OUT%FHNDL+4)
           WRITE(OUT%FHNDL+5)  SNGL(TIME)
           WRITE(OUT%FHNDL+5)  (SNGL(CURR_GLOBAL(IP,1)), SNGL(CURR_GLOBAL(IP,2)), SNGL(CURR_GLOBAL(IP,5))  , IP = 1, NP_GLOBAL)
           FLUSH(OUT%FHNDL+5)
           WRITE(OUT%FHNDL+6)  SNGL(TIME)
           WRITE(OUT%FHNDL+6)  (SNGL(WIND_GLOBAL(IP,9)), SNGL(WIND_GLOBAL(IP,8)), SNGL(WIND_GLOBAL(IP,7))  , IP = 1, NP_GLOBAL)
           FLUSH(OUT%FHNDL+6)
           WRITE(OUT%FHNDL+7)  SNGL(TIME) 
           WRITE(OUT%FHNDL+7)  (SNGL(WIND_GLOBAL(IP,4)), SNGL(WIND_GLOBAL(IP,5)), SNGL(WIND_GLOBAL(IP,6))  , IP = 1, NP_GLOBAL)
           FLUSH(OUT%FHNDL+7)
           IF (LQSTEA .AND. LCHKCONV) THEN
             WRITE(OUT%FHNDL+8) SNGL(TIME) 
             WRITE(OUT%FHNDL+8)  (SNGL(ITER_GLOBAL(IP)), SNGL(ITER_GLOBAL(IP)), SNGL(ITER_GLOBAL(IP))  , IP = 1, NP_GLOBAL)
             FLUSH(OUT%FHNDL+8)
           ENDIF
           IF (DoAirSea) THEN
             DO IP = 1, NP_GLOBAL
               WRITE(OUT%FHNDL+9,'(10F15.6)') SNGL(WIND_GLOBAL(IP,:))
             ENDDO
             FLUSH(OUT%FHNDL+9)
           END IF
           IF (LCFL) THEN
             WRITE(OUT%FHNDL+10)  SNGL(TIME)
             WRITE(OUT%FHNDL+10)  (SNGL(CFLCXY(1,IP)), SNGL(CFLCXY(2,IP)), SNGL(CFLCXY(3,IP)), IP = 1, MNP)
           ENDIF
         END IF ! myrank
#else

!$OMP DO PRIVATE (IP,ACLOC,OUTPARS,CURRPARS,WINDPARS)
         DO IP = 1, MNP
            OUTT(IP,:) = 0.
            IF (DEP(IP) .GT. DMIN) THEN
              ACLOC(:,:) = AC2(:,:,IP)
              CALL INTPAR(IP, MSC, ACLOC, OUTPARS)
              CALL CURRPAR(IP, CURRPARS)
              CALL WINDPAR(IP, WINDPARS)
            ELSE
              OUTPARS  = 0.
              CURRPARS = 0.
              WINDPARS = 0.
            END IF
            OUTT(IP,:) = OUTPARS(:)
            CURR(IP,:) = CURRPARS(:)
            WIND(IP,:) = WINDPARS(:)
         END DO

         IF (LMONO_OUT) THEN
!$OMP WORKSHARE
           OUTT(:,1) = OUTT(:,1) / SQRT(2.)
!$OMP END WORKSHARE
         ENDIF

!$OMP MASTER
         IF (LINIT_OUTPUT) THEN
           OPEN(OUT%FHNDL+1, FILE  = 'ergzusw.bin'  , FORM = 'UNFORMATTED')
           OPEN(OUT%FHNDL+2, FILE  = 'erguvh.bin'   , FORM = 'UNFORMATTED')
           OPEN(OUT%FHNDL+3, FILE  = 'ergwind.bin'  , FORM = 'UNFORMATTED')
           OPEN(OUT%FHNDL+4, FILE  = 'ergchar.bin'  , FORM = 'UNFORMATTED')
           OPEN(OUT%FHNDL+5, FILE  = 'ergtm02.bin'  , FORM = 'UNFORMATTED')
           OPEN(OUT%FHNDL+6, FILE  = 'ergshallow.bin'  , FORM = 'UNFORMATTED')
           OPEN(OUT%FHNDL+7, FILE  = 'ergufric.bin'  , FORM = 'UNFORMATTED')
           OPEN(OUT%FHNDL+8, FILE  = 'ergtau.bin'  , FORM = 'UNFORMATTED')
           IF (DoAirSea) THEN
             OPEN(OUT%FHNDL+9, FILE  = 'airsea.dat'  , FORM = 'FORMATTED')
           END IF
           OPEN(OUT%FHNDL+10, FILE  = 'cflcxy.bin'  , FORM = 'UNFORMATTED')
         END IF
         WRITE(OUT%FHNDL+1) SNGL(TIME) 
         WRITE(OUT%FHNDL+1)  (SNGL(OUTT(IP,7)), SNGL(OUTT(IP,8)), SNGL(OUTT(IP,1)), IP = 1, MNP)
         FLUSH(OUT%FHNDL+1)
         WRITE(OUT%FHNDL+2) SNGL(TIME)
         WRITE(OUT%FHNDL+2)  (SNGL(CURR(IP,1)), SNGL(CURR(IP,2)), SNGL(DEP(IP)), IP = 1, MNP)
         FLUSH(OUT%FHNDL+2)
         WRITE(OUT%FHNDL+3) SNGL(TIME)
         WRITE(OUT%FHNDL+3)  (SNGL(OUTT(IP,7)), SNGL(OUTT(IP,8)), SNGL(WIND(IP,3)), IP = 1, MNP)
         FLUSH(OUT%FHNDL+3)
         WRITE(OUT%FHNDL+4) SNGL(TIME)
         WRITE(OUT%FHNDL+4)  (SNGL(UFRIC(IP)), SNGL(Z0(IP)), SNGL(ALPHA_CH(IP)), IP = 1, MNP)
         FLUSH(OUT%FHNDL+4)
         WRITE(OUT%FHNDL+5) SNGL(TIME)
         WRITE(OUT%FHNDL+5)  (SNGL(OUTT(IP,1)), SNGL(OUTT(IP,2)), SNGL(OUTT(IP,3)), IP = 1, MNP)
         FLUSH(OUT%FHNDL+5)
         WRITE(OUT%FHNDL+6)  SNGL(TIME)
         WRITE(OUT%FHNDL+6)  (SNGL(OUTT(IP,1)), SNGL(OUTT(IP,2)), REAL(ISHALLOW(IP)), IP = 1, MNP)
         FLUSH(OUT%FHNDL+6)
         WRITE(OUT%FHNDL+7)  SNGL(TIME)
         WRITE(OUT%FHNDL+7)  (SNGL(WIND(IP,8)), SNGL(WIND(IP,9)), SNGL(WIND(IP,8)), IP = 1, MNP)
         FLUSH(OUT%FHNDL+7)
         IF (LQSTEA .AND. LCHKCONV) THEN
           WRITE(OUT%FHNDL+8)  SNGL(TIME) 
           WRITE(OUT%FHNDL+8)  (REAL(IP_IS_STEADY(IP)), REAL(IP_IS_STEADY(IP)), REAL(IP_IS_STEADY(IP))  , IP = 1, NP_TOTAL)
           FLUSH(OUT%FHNDL+8)
         ENDIF
         IF (DoAirSea) THEN
           DO IP = 1, MNP
             WRITE(OUT%FHNDL+9,'(10F15.6)') SNGL(WIND(IP,:))
           ENDDO
           FLUSH(OUT%FHNDL+9)
         END IF
         IF (LCFL) THEN
           WRITE(OUT%FHNDL+10)  SNGL(TIME)
           WRITE(OUT%FHNDL+10)  (SNGL(CFLCXY(1,IP)), SNGL(CFLCXY(2,IP)), SNGL(CFLCXY(3,IP)), IP = 1, MNP)
         ENDIF
         
!$OMP END MASTER
#endif
        WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH XFN_HISTORY'
        FLUSH(STAT%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PRINT_HS_TRIPLE(eX, eY)
      USE DATAPOOL
      implicit none
      REAL(rkind), intent(in) :: eX, eY
      REAL(rkind) :: HS, TM01, TM02, TM10, KLM, WLM
      REAL(rkind) :: X(3), Y(3), WI(3)
      REAL(rkind) :: HSinterp, HSinterpB, sumAC
      integer :: IP, I, ISMAX, IEfind
      REAL(rkind) :: ACLOC(MSC,MDC)
      REAL(rkind) :: XYTMP(2,MNP)
      XYTMP(1,:) = XP
      XYTMP(2,:) = YP
      IEfind=0
      CALL FIND_ELE(MNE,MNP,INE,XYTMP,eX,eY,IEfind)
      DO I=1,3
        IP=INE(I,IEfind)
        X(I)=XP(IP)
        Y(I)=YP(IP)
      END DO
      CALL INTELEMENT_COEF(X,Y,eX,eY,WI)
      HSinterp=0
      HSinterpB=0
      ISMAX=MSC
      sumAC=0
      DO I=1,3
        IP=INE(I,IEfind)
        ACLOC(:,:) = AC2(:,:,IP)
        CALL MEAN_PARAMETER(IP,ACLOC,ISMAX,HS,TM01,TM02,TM10,KLM,WLM)
        WRITE(DBG%FHNDL,*) 'IP=', IP, 'HS=', HS, 'wi=', WI(I)
        HSinterp=HSinterp+WI(I)*HS
        HSinterpB=HSinterpB + WI(I)*HS*HS
        sumAC=sumAC + WI(I)*sum(ACLOC)
      END DO
      HSinterpB=SQRT(HSinterpB)
      WRITE(DBG%FHNDL,*) 'HSinterp=', HSinterp, ' HS(b)=', HSinterpB
      WRITE(DBG%FHNDL,*) 'sumAC=', sumAC
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OUTPUT_STE(CTIME,LINIT_OUTPUT)
      USE DATAPOOL
      IMPLICIT NONE

      CHARACTER(LEN=15), INTENT(IN) :: CTIME
      LOGICAL, INTENT(IN)           :: LINIT_OUTPUT

      REAL(rkind) :: ACLOC(MSC,MDC), ACOUT_1D(MSC,3), ACOUT_2D(MSC*MDC)

      CHARACTER(LEN=20) :: TITLEFORMAT,OUTPUTFORMAT
      CHARACTER(LEN=2)  :: CHRTMP
      character(len=256) :: FILEWRITE
      INTEGER           :: I, IP, NI(3), IS
      LOGICAL           :: ALIVE, LSAME
      REAL(rkind)       :: WI(3)
#ifdef MPI_PARALL_GRID
      REAL(rkind) :: TheIsumR
#endif
      REAL(rkind) :: DEPLOC, WKLOC(MSC), CURTXYLOC(2)
#ifndef MPI_PARALL_GRID
      REAL(rkind) :: WATLOC, ESUM
#endif
#ifndef MPI_PARALL_GRID
      REAL(rkind) :: USTARLOC, Z0LOC, ALPHALOC, WINDXLOC, WINDYLOC, CDLOC
#endif

      WRITE(CHRTMP,'(I2)') OUTVARS
      TITLEFORMAT  = '(A15,X,X,'//TRIM(CHRTMP)//'A15)'
      OUTPUTFORMAT = '(A15,X,X,'//TRIM(CHRTMP)//'F15.8)'

      LSAME = .FALSE.

#ifndef MPI_PARALL_GRID
      DO I = 1, IOUTS ! Loop over stations ...
        IF (STATION(I)%IFOUND == 1) THEN
          STATION(I)%OUTPAR_NODE = 0.
          CALL INTELEMENT_AC_LOC(I,ACLOC,CURTXYLOC,DEPLOC,WATLOC,WKLOC)
          CALL INTPAR_LOC(I, STATION(I)%ISMAX,WKLOC,DEPLOC,CURTXYLOC,ACLOC,STATION(I)%OUTPAR_NODE)
!          WRITE(DBG%FHNDL,*) 'HS=', STATION(I)%OUTPAR_NODE(1), 'sumAC=', sum(ACLOC)
!          CALL PRINT_HS_TRIPLE(STATION(I)%XCOORD, STATION(I)%YCOORD)
          NI = INE(:,STATION(I)%ELEMENT)
          CALL INTELEMENT(XP(NI),YP(NI),UFRIC(NI), STATION(I)%XCOORD,STATION(I)%YCOORD,WI,USTARLOC,LSAME)
          CALL INTELEMENT(XP(NI),YP(NI),Z0(NI),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,Z0LOC,LSAME)
          CALL INTELEMENT(XP(NI),YP(NI),ALPHA_CH(NI),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,ALPHALOC,LSAME)
          CALL INTELEMENT(XP(NI),YP(NI),WINDXY(NI,1),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,WINDXLOC,LSAME)
          CALL INTELEMENT(XP(NI),YP(NI),WINDXY(NI,2),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,WINDYLOC,LSAME)
          CALL INTELEMENT(XP(NI),YP(NI),CD(NI),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,CDLOC,LSAME)
        ELSE
          USTARLOC = 0.
          Z0LOC = 0.
          ALPHALOC = 0.
          WINDXLOC = 0.
          WINDYLOC = 0.
          CDLOC = 0.
          WKLOC = 0.
          DEPLOC = 0.
          CURTXYLOC = 0.
          ACLOC = 0.
          WATLOC = 0.
        END IF
        ACLOC_STATIONS(:,:,I)=ACLOC
        CDLOC_STATIONS(I)=CDLOC
        Z0LOC_STATIONS(I)=Z0LOC
        ALPHALOC_STATIONS(I)=ALPHALOC
        WINDXLOC_STATIONS(I)=WINDXLOC
        WINDYLOC_STATIONS(I)=WINDYLOC
        USTARLOC_STATIONS(I)=USTARLOC
        DEPLOC_STATIONS(I)=DEPLOC
        WKLOC_STATIONS(I,:)=WKLOC
        CURTXYLOC_STATIONS(I,:)=CURTXYLOC
        WATLEVLOC_STATIONS(I)=WATLOC

        STATION(I)%OUTPAR_NODE(26) = USTARLOC_STATIONS(I)
        STATION(I)%OUTPAR_NODE(27) = Z0LOC_STATIONS(I)
        STATION(I)%OUTPAR_NODE(28) = ALPHALOC_STATIONS(I)
        STATION(I)%OUTPAR_NODE(29) = WINDXLOC_STATIONS(I)
        STATION(I)%OUTPAR_NODE(30) = WINDYLOC_STATIONS(I)
        STATION(I)%OUTPAR_NODE(31) = CDLOC_STATIONS(I)
        STATION(I)%OUTPAR_NODE(32) = CURTXYLOC_STATIONS(I,1)
        STATION(I)%OUTPAR_NODE(33) = CURTXYLOC_STATIONS(I,2)
        STATION(I)%OUTPAR_NODE(34) = DEPLOC_STATIONS(I)
        STATION(I)%OUTPAR_NODE(35) = WATLEVLOC_STATIONS(I)
      END DO
!
#else
!
      DO I = 1, IOUTS ! Loop over stations ... all threads 
        IF (STATION(I)%IFOUND .GT. 0) THEN
          CALL INTELEMENT_AC_LOC(I,ACLOC_STATIONS(:,:,I),CURTXYLOC_STATIONS(I,:),DEPLOC_STATIONS(I),WATLEVLOC_STATIONS(I),WKLOC_STATIONS(I,:))
          NI = INE(:,STATION(I)%ELEMENT)
          CALL INTELEMENT(XP(NI),YP(NI),UFRIC(NI),STATION(I)%XCOORD, STATION(I)%YCOORD,WI,USTARLOC_STATIONS(I),LSAME) 
          CALL INTELEMENT(XP(NI),YP(NI),Z0(NI),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,Z0LOC_STATIONS(I),LSAME)
          CALL INTELEMENT(XP(NI),YP(NI),ALPHA_CH(NI),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,ALPHALOC_STATIONS(I),LSAME)
          CALL INTELEMENT(XP(NI),YP(NI),WINDXY(NI,1),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,WINDXLOC_STATIONS(I),LSAME)
          CALL INTELEMENT(XP(NI),YP(NI),WINDXY(NI,2),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,WINDYLOC_STATIONS(I),LSAME)
          CALL INTELEMENT(XP(NI),YP(NI),CD(NI),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,CDLOC_STATIONS(I),LSAME)
        ELSE
          ACLOC_STATIONS(:,:,I) = 0.
          DEPLOC_STATIONS(I) = 0.
          WKLOC_STATIONS(I,:) = 0.
          USTARLOC_STATIONS(I) = 0.
          Z0LOC_STATIONS(I) = 0.
          ALPHALOC_STATIONS(I) = 0.
          WINDXLOC_STATIONS(I) = 0.
          WINDYLOC_STATIONS(I) = 0.
          CDLOC_STATIONS(I) = 0.
          WATLEVLOC_STATIONS(I) = 0.
          CURTXYLOC_STATIONS(I,:) = 0.
        END IF
      END DO

!         WRITE(STAT%FHNDL,*) 'DEPTH OF THE FOUND STATIONS', DEPLOC_STATIONS
!         FLUSH(DBG%FHNDL)

      DEPLOC_SUM = 0.
      WATLEVLOC_SUM = 0.
      CURTXYLOC_SUM = 0.
      USTAR_SUM = 0.
      Z0_SUM = 0.
      ALPHA_SUM = 0.
      CD_SUM = 0.
      WINDX_SUM = 0.
      WINDY_SUM = 0.
      WKLOC_SUM = 0.
      ACLOC_SUM = 0.

      DO I = 1, IOUTS
        CALL MPI_REDUCE(   DEPLOC_STATIONS(I),DEPLOC_SUM(I),1,rtype,MPI_SUM,0,COMM,IERR)
        CALL MPI_REDUCE(   WATLEVLOC_STATIONS(I),WATLEVLOC_SUM(I),1,rtype,MPI_SUM,0,COMM,IERR)
        CALL MPI_REDUCE(CURTXYLOC_STATIONS(I,1),CURTXYLOC_SUM(I,1),1,rtype,MPI_SUM,0,COMM,IERR)
        CALL MPI_REDUCE(CURTXYLOC_STATIONS(I,2),CURTXYLOC_SUM(I,2),1,rtype,MPI_SUM,0,COMM,IERR)
        CALL MPI_REDUCE(   USTARLOC_STATIONS(I),USTAR_SUM(I),1,rtype,MPI_SUM,0,COMM,IERR)
        CALL MPI_REDUCE(   Z0LOC_STATIONS(I),Z0_SUM(I),1,rtype,MPI_SUM,0,COMM,IERR)
        CALL MPI_REDUCE(   ALPHALOC_STATIONS(I),ALPHA_SUM(I),1,rtype,MPI_SUM,0,COMM,IERR)
        CALL MPI_REDUCE(   CDLOC_STATIONS(I),CD_SUM(I),1,rtype,MPI_SUM,0,COMM,IERR)
        CALL MPI_REDUCE(WINDXLOC_STATIONS(I),WINDX_SUM(I),1,rtype,MPI_SUM,0,COMM,IERR)
        CALL MPI_REDUCE(WINDYLOC_STATIONS(I),WINDY_SUM(I),1,rtype,MPI_SUM,0,COMM,IERR)
        CALL MPI_REDUCE(WKLOC_STATIONS(I,:),   WKLOC_SUM(I,:),MSC,rtype,MPI_SUM,0,COMM,IERR)
        CALL MPI_REDUCE(ACLOC_STATIONS(:,:,I),ACLOC_SUM(:,:,I),MSC*MDC,rtype,MPI_SUM,0,COMM,IERR)
      END DO

      IF (MYRANK == 0) THEN
        DO I = 1, IOUTS
          IF (STATION(I)%ISUM .EQ. 0) THEN
            DEPLOC_STATIONS(I)           = -999.
            CURTXYLOC_STATIONS(I,:)      = -999.
            STATION(I)%OUTPAR_NODE(1:OUTVARS) = -999.
            ACLOC_STATIONS(:,:,I)        = -999.
            WKLOC_STATIONS(I,:)          = -999.
            USTARLOC_STATIONS(I)         = -999.
            ALPHALOC_STATIONS(I)         = -999.
            Z0LOC_STATIONS(I)            = -999.
            CDLOC_STATIONS(I)            = -999.
            WINDXLOC_STATIONS(I)         = -999.
            WINDYLOC_STATIONS(I)         = -999.
            WRITE(DBG%FHNDL,*) 'STATION OUT OF MESH', I
          ELSE
            TheIsumR=MyREAL(STATION(I)%ISUM)
            DEPLOC_STATIONS(I)       = DEPLOC_SUM(I)      / TheIsumR
            CURTXYLOC_STATIONS(I,:)  = CURTXYLOC_SUM(I,:) / TheIsumR
            CURTXYLOC_STATIONS(I,:)  = CURTXYLOC_SUM(I,:) / TheIsumR
            WATLEVLOC_STATIONS(I)    = WATLEVLOC_SUM(I)   / TheIsumR
            WKLOC_STATIONS(I,:)      = WKLOC_SUM(I,:)     / TheIsumR
            ACLOC_STATIONS(:,:,I)    = ACLOC_SUM(:,:,I)   / TheIsumR
            USTARLOC_STATIONS(I)     = USTAR_SUM(I)       / TheIsumR
            Z0LOC_STATIONS(I)        = Z0_SUM(I)          / TheIsumR
            ALPHALOC_STATIONS(I)     = USTAR_SUM(I)       / TheIsumR
            CDLOC_STATIONS(I)        = CD_SUM(I)          / TheIsumR
            WINDXLOC_STATIONS(I)     = WINDX_SUM(I)       / TheIsumR
            WINDYLOC_STATIONS(I)     = WINDY_SUM(I)       / TheIsumR
            !WRITE(STAT%FHNDL,*) 'SUM BEFORE INTPAR', SUM(ACLOC_STATIONS(:,:,I))
            CALL INTPAR_LOC(I ,STATION(I)%ISMAX,WKLOC_STATIONS(I,:), DEPLOC_STATIONS(I), CURTXYLOC_STATIONS(I,:), ACLOC_STATIONS(:,:,I), STATION(I)%OUTPAR_NODE)
            STATION(I)%OUTPAR_NODE(26) = USTARLOC_STATIONS(I)
            STATION(I)%OUTPAR_NODE(27) = Z0LOC_STATIONS(I)
            STATION(I)%OUTPAR_NODE(28) = ALPHALOC_STATIONS(I)
            STATION(I)%OUTPAR_NODE(29) = WINDXLOC_STATIONS(I)
            STATION(I)%OUTPAR_NODE(30) = WINDYLOC_STATIONS(I)
            STATION(I)%OUTPAR_NODE(31) = CDLOC_STATIONS(I)
            STATION(I)%OUTPAR_NODE(32) = CURTXYLOC_STATIONS(I,1)
            STATION(I)%OUTPAR_NODE(33) = CURTXYLOC_STATIONS(I,2)
            STATION(I)%OUTPAR_NODE(34) = DEPLOC_STATIONS(I)
            STATION(I)%OUTPAR_NODE(35) = WATLEVLOC_STATIONS(I)
          END IF
        END DO
      END IF
#endif
#ifdef MPI_PARALL_GRID
      IF (MYRANK == 0) THEN
#endif
        DO I = 1, IOUTS
          ACLOC=ACLOC_STATIONS(:,:,I)
!          WRITE(STAT%FHNDL,*) I, 'SUM ACLOC 1', SUM(ACLOC)
          WKLOC=WKLOC_STATIONS(I,:)
          DEPLOC=DEPLOC_STATIONS(I)
          CURTXYLOC=CURTXYLOC_STATIONS(I,:)
          FILEWRITE=TRIM(STATION(I)%NAME)//'.site'
          INQUIRE(FILE=FILEWRITE,EXIST=ALIVE)
          IF (LINIT_OUTPUT) THEN
            OPEN(OUTPARM%FHNDL,FILE=FILEWRITE, STATUS = 'UNKNOWN')
            WRITE(OUTPARM%FHNDL, TITLEFORMAT) 'TIME', OUTT_VARNAMES
          ELSE
            IF (ALIVE.eqv..false.) THEN
              CALL WWM_ABORT('Unexisting .site file please debug')
            ENDIF
            OPEN(OUTPARM%FHNDL,FILE=FILEWRITE, STATUS = 'OLD' , POSITION = 'APPEND')
          ENDIF
          WRITE(OUTPARM%FHNDL,OUTPUTFORMAT) CTIME, STATION(I)%OUTPAR_NODE(1:OUTVARS)
          FLUSH(OUTPARM%FHNDL)
          CLOSE(OUTPARM%FHNDL)

          IF (LSP1D .OR. LSP2D) THEN
            CALL CLSPEC( WKLOC, DEPLOC, CURTXYLOC, ACLOC, ACOUT_1D, ACOUT_2D )
          END IF

          IF (LSP1D) THEN
            FILEWRITE=TRIM(STATION(I)%NAME)//'.sp1d'
            INQUIRE(FILE=FILEWRITE,EXIST=ALIVE)
            IF (LINIT_OUTPUT) THEN
              OPEN(OUTSP1D%FHNDL,FILE=FILEWRITE, STATUS = 'UNKNOWN')
              WRITE(OUTSP1D%FHNDL,*) MSC
              WRITE(OUTSP1D%FHNDL,*) MDC
              WRITE(OUTSP1D%FHNDL,*) SPSIG
              WRITE(OUTSP1D%FHNDL,*) SPDIR
#ifdef MPI_PARALL_GRID
              WRITE(OUTSP1D%FHNDL,*) STATION(I)%ISUM
#else
              WRITE(OUTSP1D%FHNDL,*) STATION(I)%IFOUND
#endif
            ELSE
              IF (ALIVE.eqv..false.) THEN
                CALL WWM_ABORT('Unexisting .sp1d file please debug')
              ENDIF
              OPEN(OUTSP1D%FHNDL,FILE=FILEWRITE, STATUS = 'OLD' , POSITION = 'APPEND')
            END IF
            WRITE(OUTSP1D%FHNDL,*) CTIME
            WRITE(OUTSP1D%FHNDL,*) DEPLOC
            WRITE(OUTSP1D%FHNDL,*) CURTXYLOC
            DO IS = 1, MSC
              WRITE(OUTSP1D%FHNDL,'(F15.8,3F20.10)') SPSIG(IS)/PI2,  ACOUT_1D(IS,1), ACOUT_1D(IS,2), ACOUT_1D(IS,3)
            END DO
            FLUSH(OUTSP1D%FHNDL)
            CLOSE(OUTSP1D%FHNDL)
          END IF

          IF (LSP2D) THEN
            FILEWRITE=TRIM(STATION(I)%NAME)//'.sp2d'
            INQUIRE(FILE=FILEWRITE,EXIST=ALIVE )
            IF (LINIT_OUTPUT) THEN
              OPEN(OUTSP2D%FHNDL,FILE=FILEWRITE, STATUS = 'UNKNOWN', FORM = 'UNFORMATTED')
              WRITE(OUTSP2D%FHNDL) MSC
              WRITE(OUTSP2D%FHNDL) MDC
              WRITE(OUTSP2D%FHNDL) SPSIG
              WRITE(OUTSP2D%FHNDL) SPDIR
#ifdef MPI_PARALL_GRID
              WRITE(OUTSP2D%FHNDL) STATION(I)%ISUM
#else
              WRITE(OUTSP2D%FHNDL) STATION(I)%IFOUND
#endif
            ELSE
              IF (ALIVE.eqv..false.) THEN
                CALL WWM_ABORT('Unexisting .sp2d file please debug')
              ENDIF
              OPEN(OUTSP2D%FHNDL,FILE=FILEWRITE, STATUS = 'OLD' , POSITION = 'APPEND', FORM = 'UNFORMATTED')
            END IF
            WRITE(OUTSP2D%FHNDL) CTIME
            WRITE(OUTSP2D%FHNDL) DEPLOC
            WRITE(OUTSP2D%FHNDL) CURTXYLOC
            WRITE(OUTSP2D%FHNDL) ACLOC
            WRITE(OUTSP2D%FHNDL) ACOUT_2D
            FLUSH(OUTSP2D%FHNDL)
            CLOSE(OUTSP2D%FHNDL)
          END IF ! LSP2D
        END DO
#ifdef MPI_PARALL_GRID
      END IF
#endif
      WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH OUTPUT_STE'
      FLUSH(STAT%FHNDL)
      END SUBROUTINE
!**********************************************************************
!* The netcdf output outs the most variables and is the most          *
!* parametrizable. In order to add some new variables in output, one  *
!* needs to:                                                          *
!* ---update OUTVAR_COMPLETE (wwm_datapl.F90)                         *
!* ---update the reading of output parameter in                       *
!*    READ_HISTORY_STATION_NAMELIST (wwm_input.F90)                   *
!* ---add full name and unit in NAMEVARIABLE (wwm_netcdf.F90)         *
!* ---add needed corrections in DETERMINE_NEEDED_COMPUTATION          *
!*      (wwm_netcdf.F90)                                              *
!* ---add the commputation of the output variable in                  *
!*    PAR_COMPLETE (history)    and PAR_COMPLETE_LOC (station)        *
!* ---wwminput.nmlref for the manual.                                 *
!**********************************************************************
#ifdef NCDF
      SUBROUTINE OUTPUT_STATION_NC
      USE NETCDF
      USE DATAPOOL
      IMPLICIT NONE

# ifdef MPI_PARALL_GRID
      REAL(rkind) :: OUTPAR_STATIONS_SUM(IOUTS,OUTVARS_COMPLETE)
      REAL(rkind) :: WK_STATIONS_SUM(IOUTS,MSC)
      REAL(rkind) :: AC_STATIONS_SUM(IOUTS,MSC,MDC)
      REAL(rkind), allocatable  :: ACOUT_1D_STATIONS_SUM(:,:,:)
      REAL(rkind), allocatable  :: ACOUT_2D_STATIONS_SUM(:,:,:)
# endif
      REAL(rkind) :: OUTPAR_STATIONS(IOUTS,OUTVARS_COMPLETE)
      REAL(rkind) :: WK_STATIONS(IOUTS,MSC)
      REAL(rkind) :: AC_STATIONS(IOUTS,MSC,MDC)
      REAL(rkind), allocatable  :: ACOUT_1D_STATIONS(:,:,:)
      REAL(rkind), allocatable  :: ACOUT_2D_STATIONS(:,:,:)
      REAL(rkind) :: eTimeDay
      REAL(rkind) :: ACLOC(MSC,MDC)
      REAL(rkind) :: ACOUT_1D(MSC,3)
      REAL(rkind) :: ACOUT_2D(MSC,MDC)
      INTEGER     :: LPOS
      character(len =256) :: FILE_NAME, PRE_FILE_NAME
      integer :: iret, ncid
      integer :: var_id
      integer :: I, irec_dim, IELOC, ISMAX
      character (len = *), parameter :: CallFct = "OUTPUT_STATION_NC"
      character (len = *), parameter :: UNITS = "units"
      character (len = *), parameter :: FULLNAME = "full-name"
      integer, save ::  recs_stat
      integer, save ::  recs_stat2
      integer, save ::  ifile = 1
      logical, save :: IsInitDone = .FALSE.
      REAL(rkind)  :: OUTPAR(OUTVARS_COMPLETE)
      character(len=40) :: eStr, eStrUnit
      character(len=80) :: eStrFullName
      REAL(rkind)       :: WI(3)
      integer eInt(1)
      integer POSITION_BEFORE_POINT
# ifdef MPI_PARALL_GRID
      REAL(rkind) :: TheIsumR
# endif
      REAL(rkind) :: DEPLOC, WATLOC, WKLOC(MSC), CURTXYLOC(2)

      IF (VAROUT_STATION%ACOUT_1D.or.VAROUT_STATION%ACOUT_2D) THEN
        allocate(ACOUT_1D_STATIONS(IOUTS, MSC, 3), ACOUT_2D_STATIONS(IOUTS, MSC, MDC), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_output, allocate error 1')
      ENDIF
      DO I = 1, IOUTS
        IF (STATION(I)%IFOUND == 1) THEN
          STATION(I)%OUTPAR_NODE = 0.
          CALL INTELEMENT_AC_LOC(I, ACLOC,CURTXYLOC, DEPLOC,WATLOC,WKLOC)
          IELOC=STATION(I)%ELEMENT
          ISMAX=STATION(I)%ISMAX
          WI=STATION(I)%WI(3)
          CALL PAR_COMPLETE_LOC(ISMAX, IELOC, WI,WKLOC, DEPLOC, CURTXYLOC, ACLOC, OUTPAR)
          IF (VAROUT_STATION%ACOUT_1D.or.VAROUT_STATION%ACOUT_2D) THEN
            CALL CLSPEC(WKLOC,DEPLOC,CURTXYLOC,ACLOC,ACOUT_1D,ACOUT_2D)
          END IF
        ELSE
          WKLOC = 0.
          ACLOC = 0.
          OUTPAR = 0.
          ACOUT_1D=0.
          ACOUT_2D=0.
        END IF
        OUTPAR_STATIONS(I,:) = OUTPAR
        WK_STATIONS(I,:) = WKLOC
        AC_STATIONS(I,:,:) = ACLOC
        IF (VAROUT_STATION%ACOUT_1D.or.VAROUT_STATION%ACOUT_2D) THEN
          ACOUT_1D_STATIONS(I,:,:)=ACOUT_1D
          ACOUT_2D_STATIONS(I,:,:)=ACOUT_2D
        END IF
      END DO
# ifdef MPI_PARALL_GRID
      IF (MULTIPLEOUT_STAT.eq.0) THEN
        OUTPAR_STATIONS_SUM=0
        WK_STATIONS_SUM=0
        AC_STATIONS_SUM=0
        IF (VAROUT_STATION%ACOUT_1D.or.VAROUT_STATION%ACOUT_2D) THEN
          allocate(ACOUT_1D_STATIONS_SUM(IOUTS,MSC,3), ACOUT_2D_STATIONS_SUM(IOUTS,MSC,MDC), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_output, allocate error 2')
          ACOUT_1D_STATIONS_SUM=0
          ACOUT_2D_STATIONS_SUM=0
        END IF
        IF (IOUTS.gt.0) THEN
          CALL MPI_REDUCE(OUTPAR_STATIONS, OUTPAR_STATIONS_SUM,         &
     &       IOUTS*OUTVARS_COMPLETE, rtype,MPI_SUM,0,comm,ierr)
          CALL MPI_REDUCE(WK_STATIONS, WK_STATIONS_SUM,                 &
     &       IOUTS*MSC, rtype,MPI_SUM,0,comm,ierr)
          CALL MPI_REDUCE(AC_STATIONS, AC_STATIONS_SUM,                 &
     &       IOUTS*MSC*MDC, rtype,MPI_SUM,0,comm,ierr)
          IF (VAROUT_STATION%ACOUT_1D.or.VAROUT_STATION%ACOUT_2D) THEN
            CALL MPI_REDUCE(ACOUT_1D_STATIONS, ACOUT_1D_STATIONS_SUM,   &
     &         IOUTS*MSC*3, rtype,MPI_SUM,0,comm,ierr)
            CALL MPI_REDUCE(ACOUT_2D_STATIONS, ACOUT_2D_STATIONS_SUM,   &
     &         IOUTS*MSC*MDC, rtype,MPI_SUM,0,comm,ierr)
          END IF
        END IF
        DO I=1,IOUTS
          TheIsumR=MyREAL(STATION(I)%ISUM)
!bug: if the stations are out of the domain we have division by zero
          OUTPAR_STATIONS(I,:)=OUTPAR_STATIONS_SUM(I,:)/TheIsumR
          WK_STATIONS(I,:) =WK_STATIONS_SUM(I,:)/TheIsumR
          AC_STATIONS(I,:,:) =AC_STATIONS_SUM(I,:,:)/TheIsumR
          IF (VAROUT_STATION%ACOUT_1D.or.VAROUT_STATION%ACOUT_2D) THEN
            ACOUT_1D_STATIONS(I,:,:)=                                   &
     &         ACOUT_1D_STATIONS_SUM(I,:,:)/TheIsumR
            ACOUT_2D_STATIONS(I,:,:)=                                   &
     &         ACOUT_2D_STATIONS_SUM(I,:,:)/TheIsumR
          ENDIF
        ENDDO
        IF (VAROUT_STATION%ACOUT_1D.or.VAROUT_STATION%ACOUT_2D) THEN
          deallocate(ACOUT_1D_STATIONS_SUM, ACOUT_2D_STATIONS_SUM)
        ENDIF
      ENDIF
# endif
      DO I=1,IOUTS
# ifdef MPI_PARALL_GRID
        IF (STATION(I)%ISUM .EQ. 0) THEN
# else
        IF (STATION(I)%IFOUND .EQ. 0) THEN
# endif
          OUTPAR_STATIONS(I,:)=-999.
          WK_STATIONS(I,:) =-999.
          AC_STATIONS(I,:,:) =-999.
          IF (VAROUT_STATION%ACOUT_1D.or.VAROUT_STATION%ACOUT_2D) THEN
            ACOUT_1D_STATIONS(I,:,:) =-999.
            ACOUT_2D_STATIONS(I,:,:) =-999.
          END IF
        END IF
      END DO
      eTimeDay=MAIN%TMJD
      LPOS=POSITION_BEFORE_POINT(OUT_STATION%FNAME)
      IF (OUT_STATION%IDEF.gt.0) THEN
         WRITE (PRE_FILE_NAME,10) OUT_STATION%FNAME(1:LPOS),ifile
  10     FORMAT (a,'_',i4.4)
      ELSE
         WRITE (PRE_FILE_NAME,20) OUT_STATION%FNAME(1:LPOS)
  20     FORMAT (a)
      ENDIF
      IF (MULTIPLEOUT_STAT.eq.0) THEN
         WRITE (FILE_NAME,30) TRIM(PRE_FILE_NAME)
  30     FORMAT (a,'.nc')
      ELSE
# ifdef MPI_PARALL_GRID
         WRITE (FILE_NAME,40) TRIM(PRE_FILE_NAME),myrank
  40     FORMAT (a,'_',i4.4,'.nc')
# endif
      ENDIF
      IF (IsInitDone.eqv..FALSE.) THEN ! At the beginning ...
        IsInitDone=.TRUE.
        recs_stat2=0
!$OMP MASTER
        IF (WriteOutputProcess_stat) THEN
          CALL DEFINE_STATION_NC(FILE_NAME, MULTIPLEOUT_STAT)
        ENDIF
!$OMP END MASTER
      END IF
      IF (LMONO_OUT) THEN
        OUTPAR_STATIONS(:,1) = OUTPAR_STATIONS(:,1) / SQRT(2.)
      ENDIF
      recs_stat2=recs_stat2 + 1
      IF (WriteOutputProcess_stat) THEN
!$OMP MASTER
        iret=nf90_open(TRIM(FILE_NAME), nf90_write, ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 29, iret)

        iret=nf90_inquire(ncid, unlimitedDimId = irec_dim)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 30, iret)

        iret=nf90_inquire_dimension(ncid, irec_dim,len = recs_stat)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 31, iret)

        recs_stat=recs_stat+1
        IF (recs_stat.ne.recs_stat2) THEN
           CALL WWM_ABORT('There are more bugs to be solved (stat)');
!bug: must be the same bug like indicated above. If the stations are not in the domain this fails ...
        ENDIF
#ifdef MPI_PARALL_GRID
        iret=nf90_inq_varid(ncid,'nproc',var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 32, iret)

        eInt(1)=nproc
        iret=nf90_put_var(ncid,var_id,eInt,start = (/1/), count=(/1/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 33, iret)
#endif
        IF (VAROUT_STATION%AC) THEN
          iret=nf90_inq_varid(ncid, 'AC', var_id)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 34, iret)

          IF (NF90_RUNTYPE == NF90_OUTTYPE_STAT) THEN
            iret=nf90_put_var(ncid,var_id,AC_STATIONS,                   &
     &        start = (/1,1,1,recs_stat/), count=(/IOUTS,MSC,MDC,1/))
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 35, iret)
          ELSE
            iret=nf90_put_var(ncid,var_id,SNGL(AC_STATIONS),             &
     &        start = (/1,1,1,recs_stat/), count=(/IOUTS,MSC,MDC,1/))
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 36, iret)
          ENDIF
        END IF
        IF (VAROUT_STATION%WK) THEN
          iret=nf90_inq_varid(ncid, 'WK', var_id)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 37, iret)
          IF (NF90_RUNTYPE == NF90_OUTTYPE_STAT) THEN
            iret=nf90_put_var(ncid,var_id,WK_STATIONS,                   &
     &        start = (/1,1,recs_stat/), count=(/IOUTS,MSC,1/))
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 38, iret)
          ELSE
            iret=nf90_put_var(ncid,var_id,SNGL(WK_STATIONS),             &
     &        start = (/1,1,recs_stat/), count=(/IOUTS,MSC,1/))
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 39, iret)
          ENDIF
        END IF
        IF (VAROUT_STATION%ACOUT_1D) THEN
          iret=nf90_inq_varid(ncid, 'ACOUT_1D', var_id)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 40, iret)
          IF (NF90_RUNTYPE == NF90_OUTTYPE_STAT) THEN
            iret=nf90_put_var(ncid,var_id,ACOUT_1D_STATIONS,             &
     &        start = (/1,1,1,recs_stat/), count=(/IOUTS,MSC,3,1/))
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 41, iret)
          ELSE
            iret=nf90_put_var(ncid,var_id,SNGL(ACOUT_1D_STATIONS),       &
     &        start = (/1,1,1,recs_stat/), count=(/IOUTS,MSC,3,1/))
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 42, iret)
          ENDIF
        END IF
        IF (VAROUT_STATION%ACOUT_2D) THEN
          iret=nf90_inq_varid(ncid, 'ACOUT_2D', var_id)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 43, iret)
          IF (NF90_RUNTYPE == NF90_OUTTYPE_STAT) THEN
            iret=nf90_put_var(ncid,var_id,ACOUT_2D_STATIONS,             &
     &        start = (/1,1,1,recs_stat/), count=(/IOUTS,MSC,MDC,1/))
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 44, iret)
          ELSE
            iret=nf90_put_var(ncid,var_id,SNGL(ACOUT_2D_STATIONS),       &
     &        start = (/1,1,1,recs_stat/), count=(/IOUTS,MSC,MDC,1/))
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 45, iret)
          ENDIF
        END IF
        DO I=1,OUTVARS_COMPLETE
          IF (VAROUT_STATION%LVAR(I)) THEN
            CALL NAMEVARIABLE(I, eStr, eStrFullName, eStrUnit)
            iret=nf90_inq_varid(ncid, TRIM(eStr), var_id)
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 46, iret)
            IF (NF90_RUNTYPE == NF90_OUTTYPE_STAT) THEN
              iret=nf90_put_var(ncid,var_id,OUTPAR_STATIONS(:,I),          &
     &              start = (/1, recs_stat/), count = (/ IOUTS, 1 /))
              CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 47, iret)
            ELSE
              iret=nf90_put_var(ncid,var_id,SNGL(OUTPAR_STATIONS(:,I)),    &
     &              start = (/1, recs_stat/), count = (/ IOUTS, 1 /))
              CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 48, iret)
            ENDIF
          END IF
        END DO
        CALL WRITE_NETCDF_TIME(ncid, recs_stat, eTimeDay)
        !write(DBG%FHNDL,*) 'Writing netcdf station record recs_stat=',recs_stat
        iret=nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 49, iret)
!$OMP END MASTER
      ENDIF
      IF (OUT_STATION%IDEF.gt.0) THEN
        IF (recs_stat2.eq.OUT_STATION%IDEF) THEN
          ifile=ifile+1
          IsInitDone = .FALSE.
        ENDIF
      ENDIF
      IF (VAROUT_STATION%ACOUT_1D.or.VAROUT_STATION%ACOUT_2D) THEN
        deallocate(ACOUT_1D_STATIONS, ACOUT_2D_STATIONS)
      ENDIF
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OUTPUT_LINE(CTIME,LINIT_OUTPUT)
         USE DATAPOOL
         IMPLICIT NONE

         CHARACTER(LEN=15), INTENT(IN) :: CTIME
         LOGICAL, INTENT(IN)           :: LINIT_OUTPUT

         REAL(rkind) :: ACLOC(MSC,MDC), ACOUT_1D(MSC,3), ACOUT_2d(MSC*MDC)
         REAL(rkind) :: TheIsumR
         CHARACTER(LEN=20) :: TITLEFORMAT,OUTPUTFORMAT
         CHARACTER(LEN=2)  :: CHRTMP
         INTEGER           :: I, IS
         LOGICAL           :: ALIVE

#ifndef MPI_PARALL_GRID
         REAL(rkind) :: DEPLOC_STATION, WATLEVLOC_STATION, WKLOC_STATION(MSC), CURTXYLOC_STATION(2), ESUM
#endif

         WRITE(CHRTMP,'(I2)') 27
         TITLEFORMAT  = '(A15,X,X,'//TRIM(CHRTMP)//'A10)'
         OUTPUTFORMAT = '(A15,X,X,'//TRIM(CHRTMP)//'F12.5)'

#ifndef MPI_PARALL_GRID
         DO I = 1, IOUTS ! Loop over stations ...

           CALL INTELEMENT_AC_LOC(I, ACLOC,CURTXYLOC_STATION, DEPLOC_STATION, WATLEVLOC_STATION,WKLOC_STATION)
           CALL INTPAR_LOC(I, STATION(I)%ISMAX,WKLOC_STATION, DEPLOC_STATION,CURTXYLOC_STATION,ACLOC,STATION(I)%OUTPAR_NODE)
           !CALL INTELEMENT_WW3GLOBAL_LOC(STATION(I)%ELEMENT,STATION(I)%XCOORD,STATION(I)%YCOORD,WW3LOCAL)

           INQUIRE(FILE=TRIM(STATION(I)%NAME)//'.site',EXIST=ALIVE)
           IF (ALIVE .AND. .NOT. LINIT_OUTPUT) THEN
             OPEN(OUTPARM%FHNDL,FILE=TRIM(STATION(I)%NAME)//'.site', STATUS = 'OLD' , POSITION = 'APPEND')
           ELSE
             OPEN(OUTPARM%FHNDL,FILE=TRIM(STATION(I)%NAME)//'.site', STATUS = 'UNKNOWN')
             WRITE(OUTPARM%FHNDL, TITLEFORMAT) 'TIME', OUTT_VARNAMES
           END IF
           WRITE(OUTPARM%FHNDL,OUTPUTFORMAT) CTIME, STATION(I)%OUTPAR_NODE
           CLOSE(OUTPARM%FHNDL)

           IF (LSP1D .OR. LSP2D) THEN
             CALL CLSPEC( WKLOC_STATION, DEPLOC_STATION, CURTXYLOC_STATION, ACLOC, ACOUT_1d, ACOUT_2d )
           END IF

           IF (LSP1D) THEN
             INQUIRE( FILE=TRIM(STATION(I)%NAME)//'.sp1d',EXIST=ALIVE)
             IF (ALIVE .AND. .NOT. LINIT_OUTPUT) THEN
               OPEN(OUTSP1D%FHNDL,FILE=TRIM(STATION(I)%NAME)//'.sp1d', STATUS = 'OLD' , POSITION = 'APPEND')
             ELSE
               OPEN(OUTSP1D%FHNDL,FILE=TRIM(STATION(I)%NAME)//'.sp1d', STATUS = 'UNKNOWN')
               WRITE(OUTSP1D%FHNDL,*) MSC, MDC
               WRITE(OUTSP1D%FHNDL,*) SPSIG, SPDIR, STATION(I)%IFOUND
             END IF
             WRITE(OUTSP1D%FHNDL,*) CTIME, WKLOC_STATION,DEPLOC_STATION, CURTXYLOC_STATION
             DO IS = 1, MSC
               WRITE(OUTSP1D%FHNDL,'(F15.8,3F20.10)') SPSIG(IS)/PI2, ACOUT_1D(IS,1), ACOUT_1D(IS,2), ACOUT_1D(IS,3)
             END DO
             CLOSE(OUTSP1D%FHNDL)
           END IF

           IF (LSP2D) THEN
             INQUIRE(FILE=TRIM(STATION(I)%NAME)//'.sp2d',EXIST=ALIVE)
             IF (ALIVE .AND. .NOT. LINIT_OUTPUT) THEN
               OPEN(OUTSP2D%FHNDL,FILE=TRIM(STATION(I)%NAME)//'.sp2d', STATUS = 'OLD' , POSITION = 'APPEND', FORM = 'UNFORMATTED')
             ELSE
               OPEN(OUTSP2D%FHNDL,FILE=TRIM(STATION(I)%NAME)//'.sp2d', STATUS = 'UNKNOWN', FORM = 'UNFORMATTED')
               WRITE(OUTSP2D%FHNDL) MSC, MDC
               WRITE(OUTSP2D%FHNDL) SPSIG, SPDIR, STATION(I)%IFOUND
             END IF
             WRITE(OUTSP2D%FHNDL) CTIME, WKLOC_STATION, DEPLOC_STATION, CURTXYLOC_STATION
             WRITE(OUTSP2D%FHNDL) ACLOC, ACOUT_2D
             CLOSE(OUTSP2D%FHNDL)
           END IF ! LSP2D

         END DO ! IOUTS
!
#else
!
         DO I = 1, IOUTS ! Loop over stations ...
           IF (STATION(I)%IFOUND .EQ. 0) CYCLE
           CALL INTELEMENT_AC_LOC(I, ACLOC_STATIONS(:,:,I), CURTXYLOC_STATIONS(I,:),DEPLOC_STATIONS(I), WATLEVLOC_STATIONS(I),WKLOC_STATIONS(I,:))
           !WRITE(DBG%FHNDL,*) 'INTERPOLATED MYRANK =', MYRANK, I, DEPLOC(I), CURTXYLOC(I,:), SUM(WKLOC(I,:)), SUM(ACLOC_STATIONS(:,:,I))
         END DO

         WRITE(DBG%FHNDL,*) 'DEPTH OF THE FOUND STATIONS LINE', DEPLOC_STATIONS

         CALL MPI_REDUCE(DEPLOC_STATIONS(:),DEPLOC_SUM(:),IOUTS,rtype,MPI_SUM,0,COMM,IERR)
         CALL MPI_REDUCE(DEPLOC_STATIONS(:),WATLEVLOC_SUM(:),IOUTS,rtype,MPI_SUM,0,COMM,IERR)
         CALL MPI_REDUCE(CURTXYLOC_STATIONS(:,1),CURTXYLOC_SUM(:,1),IOUTS,rtype,MPI_SUM,0,COMM,IERR)
         CALL MPI_REDUCE(CURTXYLOC_STATIONS(:,2),CURTXYLOC_SUM(:,2),IOUTS,rtype,MPI_SUM,0,COMM,IERR)

         DO I = 1, IOUTS
           CALL MPI_REDUCE(WKLOC_STATIONS(I,:),WKLOC_SUM(I,:),MSC,rtype,MPI_SUM,0,COMM,IERR)
           CALL MPI_REDUCE(ACLOC_STATIONS(:,:,I),ACLOC_SUM(:,:,I),MSC*MDC,rtype,MPI_SUM,0,COMM,IERR)
         END DO

         IF (MYRANK == 0) THEN

           DO I = 1, IOUTS

             INQUIRE(FILE=TRIM(STATION(I)%NAME)//'.site',EXIST=ALIVE )
             IF (ALIVE .AND. .NOT. LINIT_OUTPUT) THEN
               OPEN(OUTPARM%FHNDL,FILE=TRIM(STATION(I)%NAME)//'.site', STATUS = 'OLD' , POSITION = 'APPEND')
             ELSE
               OPEN(OUTPARM%FHNDL,FILE=TRIM(STATION(I)%NAME)//'.site', STATUS = 'UNKNOWN')
               WRITE(OUTPARM%FHNDL, TITLEFORMAT) 'TIME', OUTT_VARNAMES
             END IF

             IF (STATION(I)%ISUM .EQ. 0) THEN
               DEPLOC_STATIONS(I)       = -999.
               CURTXYLOC_STATIONS(I,:)  = -999.
               STATION(I)%OUTPAR_NODE(1:24) = -999.
               ACLOC           = -999.
               WKLOC_STATIONS           = -999.
               WRITE(DBG%FHNDL,*) 'STATION OUT OF MESH', I
             ELSE
               TheIsumR=MyREAL(STATION(I)%ISUM)
               DEPLOC_STATIONS(I)       = DEPLOC_SUM(I)      / TheIsumR
               CURTXYLOC_STATIONS(I,:)  = CURTXYLOC_SUM(I,:) / TheIsumR
               WKLOC_STATIONS(I,:)      = WKLOC_SUM(I,:)     / TheIsumR
               ACLOC           = ACLOC_SUM(:,:,I)   / TheIsumR
               CALL INTPAR_LOC(I ,STATION(I)%ISMAX,WKLOC_STATIONS(I,:), DEPLOC_STATIONS(I),CURTXYLOC_STATIONS(I,:),ACLOC, STATION(I)%OUTPAR_NODE)
             END IF

             WRITE(OUTPARM%FHNDL,OUTPUTFORMAT) CTIME, STATION(I)%OUTPAR_NODE(1:24), DEPLOC_STATIONS(I), CURTXYLOC_STATIONS(I,:)
             CLOSE(OUTPARM%FHNDL)

             IF (LSP1D .OR. LSP2D) THEN
               CALL CLSPEC( WKLOC_STATIONS(I,:), DEPLOC_STATIONS(I), CURTXYLOC_STATIONS(I,:), ACLOC, ACOUT_1d, ACOUT_2d )
             END IF

             IF (LSP1D) THEN
               INQUIRE(FILE=TRIM(STATION(I)%NAME)//'.sp1d',EXIST=ALIVE)
               IF (ALIVE .AND. .NOT. LINIT_OUTPUT) THEN
                 OPEN(OUTSP1D%FHNDL,FILE=TRIM(STATION(I)%NAME)//'.sp1d',STATUS = 'OLD' , POSITION = 'APPEND')
               ELSE
                 OPEN(OUTSP1D%FHNDL,FILE=TRIM(STATION(I)%NAME)//'.sp1d',STATUS = 'UNKNOWN')
                 WRITE(OUTSP1D%FHNDL,*) MSC, MDC
                 WRITE(OUTSP1D%FHNDL,*) SPSIG, SPDIR, STATION(I)%ISUM
               END IF
               WRITE(OUTSP1D%FHNDL,*) CTIME, WKLOC_STATIONS(I,:), DEPLOC_STATIONS(I), CURTXYLOC_STATIONS(I,:)
               DO IS = 1, MSC
                 WRITE(OUTSP1D%FHNDL,'(F15.8,3F20.10)') SPSIG(IS)/PI2, ACOUT_1D(IS,1), ACOUT_1D(IS,2), ACOUT_1D(IS,3)
               END DO
               CLOSE(OUTSP1D%FHNDL)
             END IF

             IF (LSP2D) THEN
               INQUIRE(FILE=TRIM(STATION(I)%NAME)//'.sp2d',EXIST=ALIVE)
               IF (ALIVE .AND. .NOT. LINIT_OUTPUT) THEN
                 OPEN(OUTSP2D%FHNDL,FILE=TRIM(STATION(I)%NAME)//'.sp2d', STATUS = 'OLD' , POSITION = 'APPEND', FORM = 'UNFORMATTED')
               ELSE
                 OPEN(OUTSP2D%FHNDL,FILE=TRIM(STATION(I)%NAME)//'.sp2d', STATUS = 'UNKNOWN', FORM = 'UNFORMATTED')
                 WRITE(OUTSP2D%FHNDL) MSC, MDC
                 WRITE(OUTSP2D%FHNDL) SPSIG, SPDIR, STATION(I)%ISUM
               END IF
               WRITE(OUTSP2D%FHNDL) CTIME, WKLOC_STATIONS(I,:), DEPLOC_STATIONS(I), CURTXYLOC_STATIONS(I,:)
               WRITE(OUTSP2D%FHNDL) ACLOC, ACOUT_2D
               CLOSE(OUTSP2D%FHNDL)
             END IF ! LSP2D

           END DO ! IOUTS

         END IF ! myrank
#endif
         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CURRPAR(IP, OUTPAR)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP
         REAL(rkind)   , INTENT(OUT)   :: OUTPAR(CURRVARS)

         OUTPAR    = 0.

#ifdef SCHISM
         OUTPAR(1) = UU2(NVRT,IP)         ! Current in X-direction
         OUTPAR(2) = VV2(NVRT,IP)         ! Current in Y-direction
         OUTPAR(3) = ETA2(IP)             ! Water Level
         OUTPAR(4) = ETA1(IP)             ! Water Level in last time step
         OUTPAR(5) = MAX(ZERO,WLDEP(IP) + ETA2(IP)) ! Total water depth
#else
         OUTPAR(1) = CURTXY(IP,1)         ! Current in X-direction
         OUTPAR(2) = CURTXY(IP,2)         ! Current in Y-direction
         OUTPAR(3) = WATLEV(IP)           ! Water Level
         OUTPAR(4) = WATLEVOLD(IP)        ! Water Level in last time step
         OUTPAR(5) = DEP(IP)              ! Total water depth
#endif

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WINDPAR(IP,OUTPAR)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER    , INTENT(IN)    :: IP
         REAL(rkind), INTENT(OUT)   :: OUTPAR(WINDVARS)

         OUTPAR    = 0.

         OUTPAR(1) = WINDXY(IP,1) ! wind vector u10,x
         OUTPAR(2) = WINDXY(IP,2) ! wind vector u10,y
         OUTPAR(3) = SQRT(WINDXY(IP,1)**2.+WINDXY(IP,2)**2.) ! wind magnitutde u10
         OUTPAR(4) = TAUW(IP)     ! wave stress from the discrete part of the spectra
         OUTPAR(5) = TAUHF(IP)    ! high freq. part of the waves.
         OUTPAR(6) = TAUTOT(IP)   ! total stress of the wave
         OUTPAR(7) = Z0(IP)       ! apparent rougnes lengths [m]
         OUTPAR(8) = UFRIC(IP)    ! ustar [m/s]
         OUTPAR(9) = ALPHA_CH(IP) ! Charnock Parameter gz0/ustar**2 [-}
         OUTPAR(10)= CD(IP)       ! Drag Coefficient

!         WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH WINDPAR'
!         FLUSH(STAT%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                     *
!**********************************************************************
      SUBROUTINE INTPAR(IP, ISMAX, ACLOC, OUTPAR)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)           :: IP, ISMAX
         REAL(rkind)   , INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind)   , INTENT(OUT)   :: OUTPAR(OUTVARS)

         REAL(rkind)                   :: HS,TM01,TM02,TM10,KLM,WLM
         REAL(rkind)                   :: TPP,FPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKD,DPEAK,TPPD,KPPD,CGPD,CPPD
         REAL(rkind)                   :: UBOT,ORBITAL,BOTEXPER,TMBOT,URSELL,ETOTS,ETOTC,DM,DSPR

         OUTPAR    = 0.

         CALL MEAN_PARAMETER(IP,ACLOC,ISMAX,HS,TM01,TM02,TM10,KLM,WLM)

         OUTPAR(1) = HS         ! Significant wave height
         OUTPAR(2) = TM01       ! Mean average period
         OUTPAR(3) = TM02       ! Zero down crossing period for comparison with buoy.
         OUTPAR(4) = TM10       ! Mean period of wave overtopping/run-up
         OUTPAR(5) = KLM        ! Mean wave number
         OUTPAR(6) = WLM        ! Mean wave length

!         write(DBG%FHNDL,'(6F15.8)') outpar(1:6)

         CALL MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,ISMAX,ETOTS,ETOTC,DM,DSPR)

         OUTPAR(7)  = ETOTC     ! Etot energy in horizontal direction
         OUTPAR(8)  = ETOTS     ! Etot energy in vertical direction
         OUTPAR(9)  = DM        ! Mean average energy transport direction
         OUTPAR(10) = DSPR      ! Mean directional spreading

         CALL PEAK_PARAMETER(IP,ACLOC,ISMAX,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKD,DPEAK,TPPD,KPPD,CGPD,CPPD)

         OUTPAR(11)  = TPPD     ! Discrete Peak Period
         OUTPAR(12)  = TPP      ! Peak period
         OUTPAR(13)  = CPP      ! Peak phase vel.
         OUTPAR(14)  = WNPP     ! Peak n-factor
         OUTPAR(15)  = CGPP     ! Peak group vel.
         OUTPAR(16)  = KPP      ! Peak wave number
         OUTPAR(17)  = LPP      ! Peak wave length.
         OUTPAR(18)  = PEAKD    ! Peak direction
         OUTPAR(19)  = PEAKDSPR ! Peak directional spreading
         OUTPAR(20)  = DPEAK    ! Discrete peak direction

         CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT,'INTPAR')

         OUTPAR(21)  = UBOT     ! near bottom vel. 
         OUTPAR(22)  = ORBITAL  ! Orbital vel.
         OUTPAR(23)  = BOTEXPER ! Bottom excursion period.
         OUTPAR(24)  = TMBOT    ! near bottom period. 

         IF (TPP .GT. THR) THEN
           CALL URSELL_NUMBER(HS,MyREAL(1)/TPP,DEP(IP),URSELL)
           OUTPAR(25) = URSELL  ! Ursell number based on peak period ...
         ELSE
           OUTPAR(25) = ZERO    ! Ursell number based on peak period ...
         ENDIF

         OUTPAR(26) = UFRIC(IP) ! Friction velocity 
         OUTPAR(27) = Z0(IP)    ! Roughness length
         OUTPAR(28) = ALPHA_CH(IP) ! Charnock coefficient
         OUTPAR(29) = CD(IP)    ! Drag coefficient
         OUTPAR(30) = WINDXY(IP,1) ! windx
         OUTPAR(31) = WINDXY(IP,2) ! windy

!         WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH INTPAR'
!         FLUSH(STAT%FHNDL)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEBUG_PRINT_HS_STATS
      USE DATAPOOL
      implicit none
      REAL(rkind) :: ACLOC(MSC,MDC)
      REAL(rkind) :: HS,TM01,TM02,TM10,KLM,WLM
      REAL(rkind) :: MaxHS, SumHS, AvgHS
      INTEGER IP, ISMAX
      MaxHS=ZERO
      SumHS=ZERO
      ISMAX=MSC
      DO IP=1,MNP
        ACLOC(:,:) = AC2(:,:,IP)
        CALL MEAN_PARAMETER(IP,ACLOC,ISMAX,HS,TM01,TM02,TM10,KLM,WLM)
        IF (HS.gt.MaxHS) THEN
          MaxHS=HS
        END IF
        SumHS=SumHS + HS
      END DO
      AvgHS=SumHS/MyREAL(MNP)
      WRITE(DBG%FHNDL,*) 'DEBUG AvgHS=', AvgHS, ' MaxHS=', MaxHS

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PAR_COMPLETE(IP, ISMAX, ACLOC, OUTPAR)
      USE DATAPOOL
      IMPLICIT NONE

      INTEGER, INTENT(IN)           :: IP, ISMAX
      REAL(rkind)   , INTENT(IN)    :: ACLOC(MSC,MDC)
      REAL(rkind)   , INTENT(OUT)   :: OUTPAR(OUTVARS_COMPLETE)

      REAL(rkind)                   :: HS,TM01,TM02,TM10,KLM,WLM
      REAL(rkind)                   :: TPP,FPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKD,DPEAK,TPPD,KPPD,CGPD,CPPD
      REAL(rkind)                   :: UBOT,ORBITAL,BOTEXPER,TMBOT,URSELL,ETOTS,ETOTC,DM,DSPR
      REAL(rkind)                   :: STOKESSURFX,STOKESSURFY
      REAL(rkind)                   :: STOKESBAROX,STOKESBAROY
      REAL(rkind)                   :: STOKESBOTTX,STOKESBOTTY

      OUTPAR    = 0.

      IF (VAROUT_HISTORY%ComputeMean) THEN
        CALL MEAN_PARAMETER(IP,ACLOC,ISMAX,HS,TM01,TM02,TM10,KLM,WLM)
        OUTPAR(1) = HS         ! Significant wave height
        OUTPAR(2) = TM01       ! Mean average period
        OUTPAR(3) = TM02       ! Zero down crossing period for comparison with buoy.
        OUTPAR(4) = TM10       ! Mean period of wave overtopping/run-up 
        OUTPAR(5) = KLM        ! Mean wave number
        OUTPAR(6) = WLM        ! Mean wave length
      END IF

      IF (VAROUT_HISTORY%ComputeDirSpread) THEN
        CALL MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,ISMAX,ETOTS,ETOTC,DM,DSPR)
        OUTPAR(7)   = ETOTC     ! Etot energy in horizontal direction
        OUTPAR(8)   = ETOTS     ! Etot energy in vertical direction
        OUTPAR(9)   = DM        ! Mean average energy transport direction
        OUTPAR(10)  = DSPR      ! Mean directional spreading
      END IF

      IF (VAROUT_HISTORY%ComputePeak) THEN
        CALL PEAK_PARAMETER(IP,ACLOC,ISMAX,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKD,DPEAK,TPPD,KPPD,CGPD,CPPD)
        OUTPAR(11)  = TPPD     ! Discrete Peak Period
        OUTPAR(12)  = CPPD     ! Discrete Peak speed
        OUTPAR(13)  = KPPD     ! Discrete Peak wave number
        OUTPAR(14)  = CGPD     ! Discrete Peak group speed
        OUTPAR(15)  = TPP      ! Peak period
        OUTPAR(16)  = CPP      ! Peak phase vel.
        OUTPAR(17)  = WNPP     ! Peak n-factor
        OUTPAR(18)  = CGPP     ! Peak group vel.
        OUTPAR(19)  = KPP      ! Peak wave number
        OUTPAR(20)  = LPP      ! Peak wave length.
        OUTPAR(21)  = PEAKD    ! Peak direction
        OUTPAR(22)  = PEAKDSPR ! Peak directional spreading
        OUTPAR(23)  = DPEAK    ! Discrete peak direction
      END IF

      IF (VAROUT_HISTORY%ComputeCurr) THEN
        CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT,'PAR_COMPLETE')
        OUTPAR(24)  = UBOT     !
        OUTPAR(25)  = ORBITAL  ! Orbital vel.
        OUTPAR(26)  = BOTEXPER ! Bottom excursion period.
        OUTPAR(27)  = TMBOT    !
      END IF

      IF (VAROUT_HISTORY%ComputeUrsell) THEN
        CALL URSELL_NUMBER(HS,MyREAL(1)/TPP,DEP(IP),URSELL)
        OUTPAR(28) = URSELL    ! Uresell number based on peak period ...
      END IF

      OUTPAR(29) = UFRIC(IP)
      OUTPAR(30) = Z0(IP)
      OUTPAR(31) = ALPHA_CH(IP)
      OUTPAR(32) = WINDXY(IP,1)
      OUTPAR(33) = WINDXY(IP,2)
      OUTPAR(34) = CD(IP)

      OUTPAR(35) = CURTXY(IP,1)
      OUTPAR(36) = CURTXY(IP,2)
      OUTPAR(37) = WATLEV(IP)
      OUTPAR(38) = WATLEVOLD(IP)
      OUTPAR(39) = DEPDT(IP)
      OUTPAR(40) = DEP(IP)
      OUTPAR(41) = SQRT(WINDXY(IP,1)**2.+WINDXY(IP,2)**2.)
      OUTPAR(42) = TAUW(IP)
      OUTPAR(43) = TAUW(IP)
      OUTPAR(44) = TAUW(IP)
      OUTPAR(45) = TAUHF(IP)
      OUTPAR(46) = TAUTOT(IP)

      IF (VAROUT_HISTORY%ComputeStokes) THEN
        CALL STOKES_DRIFT_SURFACE_BAROTROPIC(IP,STOKESBOTTX, STOKESBOTTY,STOKESSURFX,STOKESSURFY,STOKESBAROX,STOKESBAROY)
        OUTPAR(47) = STOKESBOTTX
        OUTPAR(48) = STOKESBOTTY
        OUTPAR(49) = STOKESSURFX
        OUTPAR(50) = STOKESSURFY
        OUTPAR(51) = STOKESBAROX
        OUTPAR(52) = STOKESBAROY
      END IF
      OUTPAR(53) = RSXX(IP)
      OUTPAR(54) = RSXY(IP)
      OUTPAR(55) = RSYY(IP)
      IF (LCFL) THEN
        OUTPAR(56) = CFLCXY(1,IP)
        OUTPAR(57) = CFLCXY(2,IP)
        OUTPAR(58) = CFLCXY(3,IP)
      ENDIF
#ifdef WWM_SETUP
      IF (LZETA_SETUP) THEN
        OUTPAR(59) = ZETA_SETUP(IP)
      END IF
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                     *
!**********************************************************************
      SUBROUTINE PAR_COMPLETE_LOC(ISMAX, IELOC, WI, WKLOC, DEPLOC, CURTXYLOC, ACLOC, OUTPAR) 
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: ISMAX
      INTEGER, intent(in)    :: IELOC
      REAL(rkind), intent(in)  :: WI(3)
      REAL(rkind), INTENT(IN)  :: ACLOC(MSC,MDC), WKLOC(MSC), DEPLOC, CURTXYLOC(2)
      REAL(rkind), INTENT(OUT) :: OUTPAR(OUTVARS_COMPLETE)
      integer :: IP, J
      REAL(rkind) :: HS,TM01,TM02,KLM,WLM,TM10
      REAL(rkind) :: TPP,FPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKD,DPEAK,TPPD,KPPD,CGPD,CPPD
      REAL(rkind) :: UBOT,ORBITAL,BOTEXPER,TMBOT,URSELL,ETOTS,ETOTC,DM,DSPR
      REAL(rkind)                   :: STOKESSURFX,STOKESSURFY
      REAL(rkind)                   :: STOKESBAROX,STOKESBAROY
      REAL(rkind)                   :: STOKESBOTTX,STOKESBOTTY
      REAL(rkind) :: eWindMag
      OUTPAR    = ZERO

      IF (VAROUT_STATION%ComputeMean) THEN
        CALL MEAN_PARAMETER_LOC(ACLOC,CURTXYLOC,DEPLOC,WKLOC,ISMAX,HS,TM01,TM02,TM10,KLM,WLM)
        OUTPAR(1) = HS         ! Significant wave height
        OUTPAR(2) = TM01       ! Mean average period
        OUTPAR(3) = TM02       ! Zero down crossing period for comparison with buoy.
        OUTPAR(4) = TM10       ! Mean period of wave overtopping/run-up 
        OUTPAR(5) = KLM        ! Mean wave number
        OUTPAR(6) = WLM        ! Mean wave length
      END IF

      IF (VAROUT_STATION%ComputeDirSpread) THEN
        CALL MEAN_DIRECTION_AND_SPREAD_LOC(ACLOC,ISMAX,ETOTS,ETOTC,DM,DSPR)
        OUTPAR(7)  = ETOTC     ! Etot energy in horizontal direction
        OUTPAR(8)  = ETOTS     ! Etot energy in vertical direction
        OUTPAR(9)  = DM        ! Mean average energy transport direction
        OUTPAR(10) = DSPR      ! Mean directional spreading
      END IF

      IF (VAROUT_STATION%ComputePeak) THEN
        CALL PEAK_PARAMETER_LOC(ACLOC,DEPLOC,ISMAX,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKD,DPEAK,TPPD,KPPD,CGPD,CPPD)
        OUTPAR(11)  = TPPD     ! Discrete Peak Period
        OUTPAR(12)  = CPPD     ! Discrete Peak speed
        OUTPAR(13)  = KPPD     ! Discrete Peak wave number
        OUTPAR(14)  = CGPD     ! Discrete Peak group speed
        OUTPAR(15)  = TPP      ! Continues Peak period
        OUTPAR(16)  = CPP      ! Peak phase vel.
        OUTPAR(17)  = WNPP     ! Peak n-factor
        OUTPAR(18)  = CGPP     ! Peak group vel.
        OUTPAR(19)  = KPP      ! Peak wave number
        OUTPAR(20)  = LPP      ! Peak wave length.
        OUTPAR(21)  = PEAKD    ! Peak direction
        OUTPAR(22)  = PEAKDSPR ! Peak directional spreading
        OUTPAR(23)  = DPEAK    ! Discrete peak direction
      END IF

      IF (VAROUT_STATION%ComputeCurr) THEN
        CALL WAVE_CURRENT_PARAMETER_LOC(ACLOC,CURTXYLOC,DEPLOC,WKLOC,UBOT,ORBITAL,BOTEXPER,TMBOT)
        OUTPAR(24)  = UBOT     !
        OUTPAR(25)  = ORBITAL  ! Orbital vel.
        OUTPAR(26)  = BOTEXPER ! Bottom excursion period.
        OUTPAR(27)  = TMBOT    !
      END IF

      IF (VAROUT_STATION%ComputeUrsell) THEN
        CALL URSELL_NUMBER(HS,MyREAL(1)/TPP,DEPLOC,URSELL)
        OUTPAR(28) = URSELL    ! Uresell number based on peak period ...
      END IF

      DO J=1,3
        IP=INE(J, IELOC)
        OUTPAR(29) = OUTPAR(29) + WI(J)*UFRIC(IP)
        OUTPAR(30) = OUTPAR(30) + WI(J)*Z0(IP)
        OUTPAR(31) = OUTPAR(31) + WI(J)*ALPHA_CH(IP)
        OUTPAR(32) = OUTPAR(22) + WI(J)*WINDXY(IP,1)
        OUTPAR(33) = OUTPAR(33) + WI(J)*WINDXY(IP,2)
        OUTPAR(34) = OUTPAR(34) + WI(J)*CD(IP)
        OUTPAR(35) = OUTPAR(35) + WI(J)*CURTXY(IP,1)
        OUTPAR(36) = OUTPAR(36) + WI(J)*CURTXY(IP,2)
        OUTPAR(37) = OUTPAR(37) + WI(J)*WATLEV(IP)
        OUTPAR(38) = OUTPAR(38) + WI(J)*WATLEVOLD(IP)
        OUTPAR(39) = OUTPAR(39) + WI(J)*DEPDT(IP)
        OUTPAR(40) = OUTPAR(40) + WI(J)*DEP(IP)
        eWindMag=SQRT(WINDXY(IP,1)**2.+WINDXY(IP,2)**2.)
        OUTPAR(41) = OUTPAR(41) + WI(J)*eWindMag
        OUTPAR(42) = OUTPAR(42) + WI(J)*TAUW(IP)
        OUTPAR(43) = OUTPAR(43) + WI(J)*TAUW(IP)
        OUTPAR(44) = OUTPAR(44) + WI(J)*TAUW(IP)
        OUTPAR(45) = OUTPAR(45) + WI(J)*TAUHF(IP)
        OUTPAR(46) = OUTPAR(46) + WI(J)*TAUTOT(IP)
        OUTPAR(53) = OUTPAR(53) + WI(J)*RSXX(IP)
        OUTPAR(54) = OUTPAR(54) + WI(J)*RSXY(IP)
        OUTPAR(55) = OUTPAR(55) + WI(J)*RSYY(IP)
        IF (LCFL) THEN
          OUTPAR(56) = OUTPAR(56) + WI(J)*CFLCXY(1,IP)
          OUTPAR(57) = OUTPAR(57) + WI(J)*CFLCXY(2,IP)
          OUTPAR(58) = OUTPAR(58) + WI(J)*CFLCXY(3,IP)
        ENDIF
#ifdef WWM_SETUP
        IF (LZETA_SETUP) THEN
          OUTPAR(59) = OUTPAR(59) + WI(J)*ZETA_SETUP(IP)
        END IF
#endif
      END DO
      IF (VAROUT_STATION%ComputeStokes) THEN
        CALL STOKES_DRIFT_SURFACE_BAROTROPIC_LOC(ACLOC,DEPLOC,WKLOC,STOKESBOTTX, STOKESBOTTY,STOKESSURFX,STOKESSURFY,STOKESBAROX,STOKESBAROY)
        OUTPAR(47) = STOKESBOTTX
        OUTPAR(48) = STOKESBOTTY
        OUTPAR(49) = STOKESSURFX
        OUTPAR(50) = STOKESSURFY
        OUTPAR(51) = STOKESBAROX
        OUTPAR(52) = STOKESBAROY
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                     *
!**********************************************************************
      SUBROUTINE INTPAR_LOC(I, ISMAX, WKLOC, DEPLOC, CURTXYLOC, ACLOC, OUTPAR)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: ISMAX, I
         REAL(rkind)   , INTENT(IN)    :: ACLOC(MSC,MDC), WKLOC(MSC), DEPLOC, CURTXYLOC(2)
         REAL(rkind)   , INTENT(OUT)   :: OUTPAR(OUTVARS)

         REAL(rkind)                   :: HS,TM01,TM02,KLM,WLM,TM10
         REAL(rkind)                   :: TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKD,DPEAK,TPPD,KPPD,CGPD,CPPD
         REAL(rkind)                   :: UBOT,ORBITAL,BOTEXPER,TMBOT,URSELL,ETOTS,ETOTC,DM,DSPR,FPP

         OUTPAR = 0.

         IF (ISMAX .GT. MSC) WRITE(DBG%FHNDL,*) 'ERROR IN ISMAX INTPAR_LOC'
         IF (DEPLOC .LT. DMIN) RETURN

         CALL MEAN_PARAMETER_LOC(ACLOC,CURTXYLOC,DEPLOC,WKLOC,ISMAX,HS,TM01,TM02,TM10,KLM,WLM)

         OUTPAR(1) = HS         ! Significant wave height
         OUTPAR(2) = TM01       ! Mean average period
         OUTPAR(3) = TM02       ! Zero down crossing period for comparison with buoy.
         OUTPAR(4) = TM10       ! Mean period of wave overtopping/run-up 
         OUTPAR(5) = KLM        ! Mean wave number
         OUTPAR(6) = WLM        ! Mean wave length

         CALL MEAN_DIRECTION_AND_SPREAD_LOC(ACLOC,ISMAX,ETOTS,ETOTC,DM,DSPR)

         OUTPAR(7)  = ETOTC     ! Etot energy in horizontal direction
         OUTPAR(8)  = ETOTS     ! Etot energy in vertical direction
         OUTPAR(9)  = DM        ! Mean average energy transport direction
         OUTPAR(10)  = DSPR      ! Mean directional spreading

         CALL PEAK_PARAMETER_LOC(ACLOC,DEPLOC,ISMAX,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKD,DPEAK,TPPD,KPPD,CGPD,CPPD)

         OUTPAR(11)  = TPPD     ! Discrete Peak Period
         OUTPAR(12)  = TPP      ! Continues Peak period
         OUTPAR(13)  = CPP      ! Peak phase vel.
         OUTPAR(14)  = WNPP     ! Peak n-factor
         OUTPAR(15)  = CGPP     ! Peak group vel.
         OUTPAR(16)  = KPP      ! Peak wave number
         OUTPAR(17)  = LPP      ! Peak wave length.
         OUTPAR(18)  = PEAKD    ! Peak direction
         OUTPAR(19)  = PEAKDSPR ! Peak directional spreading
         OUTPAR(20)  = DPEAK    ! Discrete peak direction

         CALL WAVE_CURRENT_PARAMETER_LOC(ACLOC,CURTXYLOC,DEPLOC,WKLOC,UBOT,ORBITAL,BOTEXPER,TMBOT)

         OUTPAR(21)  = UBOT     ! near bottom vel. 
         OUTPAR(22)  = ORBITAL  ! Orbital vel.
         OUTPAR(23)  = BOTEXPER ! Bottom excursion period.
         OUTPAR(24)  = TMBOT    ! near bottom period. 

         IF (TPP .GT. THR) THEN
           CALL URSELL_NUMBER(HS,MyREAL(1)/TPP,DEPLOC,URSELL)
           OUTPAR(25) = URSELL    ! Uresell number based on peak period ...
         ENDIF

         OUTPAR(26:35) = ZERO 
      END SUBROUTINE
!**********************************************************************
!*                                                                     *
!**********************************************************************
      SUBROUTINE INTSPEC(MSC, MDC, DDIR, ACLOC, SPSIG, HS)
         USE DATAPOOL, ONLY : rkind, stat
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: MSC, MDC

         REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC), SPSIG(MSC), DDIR
         REAL(rkind), INTENT(OUT)   :: HS


         INTEGER :: ID, IS
         REAL(rkind)    :: ETOT, EAD , DS

         ETOT = 0.
         DO ID = 1, MDC
            DO IS = 2, MSC
               DS = SPSIG(IS) - SPSIG(IS-1)
               EAD = 0.5*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS*DDIR
               ETOT = ETOT + EAD
            END DO
         END DO

         IF (ETOT > 0.0) THEN
           HS = 4.0*SQRT(ETOT)
         ELSE
           HS = 0.0
         END IF

!         WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH INTSPEC'
!         FLUSH(STAT%FHNDL)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OUTPUT_TEST()
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP
         REAL(rkind) :: OUTT(MNP,OUTVARS), ACLOC(MSC,MDC), OUTPAR(OUTVARS)
         CHARACTER(LEN=15) :: CTIME

         CALL MJD2CT(MAIN%TMJD, CTIME)

         OPEN(MISC%FHNDL, FILE = 'output1d.dat'  , STATUS = 'UNKNOWN' )

         DO IP = 1, MNP
            IF (DEP(IP) .GT. DMIN) THEN
              ACLOC(:,:) = AC2(:,:,IP)
              CALL INTPAR( IP, MSC, ACLOC, OUTPAR )
              OUTT(IP,:) = OUTPAR(:)
            ELSE
              OUTT(IP,:) = 0.
            END IF
         END DO

         WRITE(MISC%FHNDL,'(A12,6A14)') 'IP','XP','YP','DEPTH','HS','TM01','DSPR'

         DO IP  = 1, MNP
           IF (YP(IP) .EQ. 198000.00) THEN
           WRITE (MISC%FHNDL,'(I12,13F14.2)') IP, XP(IP), YP(IP), DEP(IP), OUTT(IP,1), OUTT(IP,3), OUTT(IP,6)
           END IF
         END DO

         CLOSE(MISC%FHNDL)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef DARKO
      SUBROUTINE OUTPUT_HISTORY_SHP( TIME )
!
!     OUTPUT DARKO CSV
!
         USE DATAPOOL
         IMPLICIT NONE
         REAL(rkind), INTENT(IN) :: TIME
         INTEGER :: IP, IE, IT
         CHARACTER(LEN=15) :: CTIME
         REAL(rkind) :: ACLOC(MSC,MDC), OUTPAR(OUTVARS)
         REAL(rkind) :: OUTT(MNP,OUTVARS)

         IF (TIME .LT. THR) THEN
           OPEN( 4001, FILE = 'depth.dat' )
           OPEN( 4002, FILE = 'hs.dat' )
         END IF

         CALL MJD2CT(MAIN%TMJD, CTIME)
         DO IP = 1, MNP
            ACLOC(:,:) = AC2(:,:,IP)
            CALL INTPAR( IP, MSC, ACLOC, OUTPAR )
            OUTT(IP,:) = OUTPAR(:)
         END DO

         IF (TIME .LT. THR) THEN
           WRITE(OUT%FHNDL,'(4A10)') 'ID', 'X', 'Y', 'Z'
           DO IE = 1, MNE
             WRITE(4001,110) IE,',',XP(INE(1,IE)),',',YP(INE(1,IE)),',',DEP(INE(1,IE))
             WRITE(4001,110) IE,',',XP(INE(2,IE)),',',YP(INE(2,IE)),',',DEP(INE(2,IE))
             WRITE(4001,110) IE,',',XP(INE(3,IE)),',',YP(INE(3,IE)),',',DEP(INE(3,IE))
           END DO
         END IF

         WRITE(4002,'(4A10)') 'ID', 'X', 'Y', 'HS'
         DO IE = 1, MNE
           WRITE(4002,110) IE,',',XP(INE(1,IE)),',',YP(INE(1,IE)),',',OUTT(INE(1,IE),1)
           WRITE(4002,110) IE,',',XP(INE(2,IE)),',',YP(INE(2,IE)),',',OUTT(INE(2,IE),1)
           WRITE(4002,110) IE,',',XP(INE(3,IE)),',',YP(INE(3,IE)),',',OUTT(INE(3,IE),1)
         END DO

110      FORMAT (2X,I10,3(A2,F15.8))

        WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH OUTPUT_HISTORY_SHP'
        FLUSH(STAT%FHNDL)
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CLSPEC( WKLOC, DEPLOC, CURTXYLOC, ACLOC, ACOUT_1D, ACOUT_2D )
         USE DATAPOOL
         IMPLICIT NONE
         REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC), WKLOC(MSC), DEPLOC, CURTXYLOC(2)
         REAL(rkind), INTENT(OUT)   :: ACOUT_1D(MSC,3), ACOUT_2D(MSC,MDC)
         INTEGER :: IS, ID, ISS

         REAL(rkind)    :: UNITFAC
         REAL(rkind)    :: UDIR, DM
         REAL(rkind)    :: ECLL, EE, RR, EADD, VEC2DEG
         REAL(rkind)    :: SIG1, C1, K1, CG1, SIG2, C2, K2, CG2
         REAL(rkind)    :: OMEG1, OMEG2, OMEGA, OMEGB, DSIG, DOMEG
         REAL(rkind)    :: RLOW, RUPP, WN1, WN2
         REAL(rkind)    :: EX, EY, FF, DEG
!
         IF (LENERGY) THEN
           UNITFAC = PWIND(2) * G9  ! for true energy
         ELSE
           UNITFAC = 1.0
         END IF

         ACOUT_1D = 0.
         ACOUT_2D = 0.
         DO ID = 1, MDC
           UDIR = CURTXYLOC(1) * COSTH(ID) +  CURTXYLOC(2) * SINTH(ID)
           DO IS = 1, MSC
             ECLL = ACLOC(IS,ID)*SPSIG(IS)
             IF (.NOT. LSTCU .AND. .NOT. LSECU) THEN ! No current ... relative freq.
               ACOUT_2D(IS,ID) = ECLL
               ACOUT_1D(IS,1) = ACOUT_1D(IS,1) + ECLL * DDIR
               ACOUT_1D(IS,2) = ACOUT_1D(IS,2) + ECLL * DDIR * COSTH(ID)
               ACOUT_1D(IS,3) = ACOUT_1D(IS,3) + ECLL * DDIR * SINTH(ID)
             ELSE                                    ! Currents absolute freq.
               SIG1 = SPSIG(IS) / FRINTH
               CALL WAVEKCG(DEPLOC, SIG1, WN1, C1, K1, CG1)
               OMEG1 = SIG1 + K1 * UDIR
               SIG2 = SPSIG(IS) * FRINTH
               CALL WAVEKCG(DEPLOC, SIG2, WN2, C2, K2, CG2)
               OMEG2 = SIG2 + K2 * UDIR
               DSIG = FRINTF * SPSIG(IS)
               EE = ECLL * DSIG / ABS(OMEG2-OMEG1) ! Jacobian
               IF (OMEG1 > OMEG2) THEN             ! Swap ...
                 RR = OMEG2
                 OMEG2 = OMEG1
                 OMEG1 = RR
               END IF
               DO ISS = 1, MSC
                 OMEGA = SPSIG(ISS) / FRINTH
                 OMEGB = SPSIG(ISS) * FRINTH
                 IF (OMEG1 < OMEGB) THEN
                   RLOW = MAX(OMEG1,OMEGA)
                 ELSE
                   CYCLE
                 END IF
                 IF (OMEG2 > OMEGA) THEN
                   RUPP = MIN(OMEG2,OMEGB)
                 ELSE
                   CYCLE
                 END IF
                 IF (RUPP < RLOW) THEN
                   CALL WWM_ABORT('ERROR IN OUTPUTSPC')
                 ELSE
                   ACOUT_2D(ISS,ID) = ACOUT_2D(ISS,ID) + EE * (RUPP-RLOW)
                   EADD = EE * DDIR * (RUPP-RLOW)
                   ACOUT_1D(ISS,1) = ACOUT_1D(ISS,1) + EADD
                   ACOUT_1D(ISS,2) = ACOUT_1D(ISS,2) + EADD * COSTH(ID)
                   ACOUT_1D(ISS,3) = ACOUT_1D(ISS,3) + EADD * SINTH(ID)
                 END IF
               END DO
             END IF
           END DO ! IS
         END DO ! ID

         IF (LSTCU .OR. LSECU) THEN
         !IF (.NOT. LSTCU .OR. .NOT. LSECU) THEN
           DO ID = 1, MDC
             DO IS = 1, MSC
               DOMEG = FRINTF * SPSIG(IS)
               ACOUT_2D(IS,ID) = ACOUT_2D(IS,ID) / DOMEG
             END DO
           END DO
           DO IS = 1, MSC
             DOMEG = FRINTF * SPSIG(IS)
             ACOUT_1D(IS,:) = ACOUT_1D(IS,:) / DOMEG
           END DO
         END IF

         DO IS = 1, MSC
           IF (ACOUT_1D(IS,1) > TINY(1.0)) THEN
             EX = ACOUT_1D(IS,2) / ACOUT_1D(IS,1)
             EY = ACOUT_1D(IS,3) / ACOUT_1D(IS,1)
             DM = VEC2DEG(EX,EY)
             CALL DEG2NAUT(DM,DEG,LNAUTIN)
! overwrite x and y energy components by deg and dspr ...
             ACOUT_1D(IS,2) = DEG
             FF = MIN(1.0_rkind,SQRT(EX**2.0+EY**2.0))
             ACOUT_1D(IS,3) = SQRT(2.0_rkind-2.0_rkind*FF)*180.0_rkind/PI
! to convert ACBIN from m^2/rad/s to m^2/Hz
             ACOUT_1D(IS,1) = ACOUT_1D(IS,1) * PI2
           ELSE
             ACOUT_1D(IS,1) = -999.0
             ACOUT_1D(IS,2) = -999.0
             ACOUT_1D(IS,3) = -999.0
           END IF
         END DO

!        WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH CLSPEC'
!        FLUSH(STAT%FHNDL)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef NCDF
      SUBROUTINE HISTORY_NC_PRINTMMA(eStr, OUTT, NPWORK, NBVAR, I)
      USE DATAPOOL, only : rkind, MULTIPLEOUT_HIS, MNP, DBG, STAT
#ifdef MPI_PARALL_GRID
      USE DATAPOOL, only : nwild_loc_res, np_global
      USE DATAPOOL, only : COMM, IERR, NPROC, myrank, rtype, istatus, ierr
#endif
      implicit none
      character(len=*) :: eStr
      integer, intent(in) :: I, NPWORK, NBVAR
      real(rkind), intent(in) :: OUTT(NPWORK, NBVAR)
      integer istat
# ifdef MPI_PARALL_GRID
      REAL(rkind), allocatable :: LVect(:,:)
      REAL(rkind) :: eVect(3)
      integer :: iProc
# endif
      REAL(rkind)  :: MinV, MaxV, AvgV
# ifdef MPI_PARALL_GRID
      IF (MULTIPLEOUT_HIS.eq.0) THEN
        MinV=minval(OUTT(:,I))
        MaxV=maxval(OUTT(:,I))
        AvgV=sum(OUTT(:,I))/np_global
      ELSE
        eVect(1)=minval(OUTT(:,I))
        eVect(2)=maxval(OUTT(:,I))
        eVect(3)=sum(OUTT(:,I)*nwild_loc_res)
        IF (myrank.eq.0) THEN
          allocate(LVect(nproc,3), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_output, allocate error 3')
          LVect(1,:)=eVect
          DO iProc=2,nproc
            CALL MPI_RECV(eVect,3,rtype, iProc-1, 190, comm, istatus, ierr)
            LVect(iProc,:)=eVect
          END DO
          MinV=minval(LVect(:,1))
          MaxV=maxval(LVect(:,2))
          AvgV=sum(LVect(:,3))/np_global
          deallocate(LVect)
        ELSE
          CALL MPI_SEND(eVect,3,rtype, 0, 190, comm, ierr)
        ENDIF
      ENDIF
      IF (myrank.eq.0) THEN
        WRITE(STAT%FHNDL,110) TRIM(eStr), MinV, MaxV, AvgV
      END IF
# else
      MinV=minval(OUTT(:,I))
      MaxV=maxval(OUTT(:,I))
      AvgV=sum(OUTT(:,I))/MNP
      WRITE(STAT%FHNDL,110) TRIM(eStr), MinV, MaxV, AvgV
# endif
110     FORMAT (a8, ' : min=', F11.5, ' max=', F11.5, ' avg=', F11.5)

!        WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH HISTORY_NC_PRINTMMA'
!        FLUSH(STAT%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OUTPUT_HISTORY_NC
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
      INTEGER            :: IP
# ifdef MPI_PARALL_GRID
      REAL(rkind), allocatable  :: OUTT_LOC(:,:)
      REAL(rkind), allocatable  :: OUTT(:,:)
# else
      REAL(rkind)               :: OUTT(MNP,OUTVARS_COMPLETE)
# endif
      REAL(rkind)               :: ACLOC(MSC,MDC), OUTPAR(OUTVARS_COMPLETE)
      REAL(rkind)  :: eTimeDay
! NC defs
      INTEGER :: LPOS
      character(len =256) :: FILE_NAME, PRE_FILE_NAME
      integer :: iret, ncid, ntime_dims, nnode_dims
      integer :: var_id
      integer :: I, irec_dim
      integer POSITION_BEFORE_POINT
      character (len = *), parameter :: CallFct="OUTPUT_HISTORY_NC"
      character (len = *), parameter :: UNITS = "units"
      character (len = *), parameter :: FULLNAME = "full-name"
      character (len = *), parameter :: COORD = "coordinates"
      character(len=1), parameter :: ePoint = '.'
      integer, save ::  recs_his
      integer, save ::  recs_his2
      integer, save ::  ifile = 1
      logical, save :: IsInitDone = .FALSE.
      integer :: np_write, ne_write
      character(len=40) :: eStr, eStrUnit
      character(len=80) :: eStrFullName
      integer, allocatable :: IOBPDoutput(:,:)
      integer nbTime
!
      eTimeDay=MAIN%TMJD
# ifdef MPI_PARALL_GRID
      IF (MULTIPLEOUT_HIS.eq.0) THEN
        np_write=np_global
        ne_write=ne_global
      ELSE
        np_write=NP_RES
        ne_write=MNE
      ENDIF
# else
      np_write=MNP
      ne_write=MNE
# endif
      IF ((np_write.eq.0).or.(ne_write.eq.0)) THEN
        CALL WWM_ABORT('OUTPUT_HISTORY_NC np_write=0 or ne_write=0')
      ENDIF
      LPOS=POSITION_BEFORE_POINT(OUT_HISTORY%FNAME)
      IF (OUT_HISTORY%IDEF.gt.0) THEN
         WRITE (PRE_FILE_NAME,10) OUT_HISTORY%FNAME(1:LPOS),ifile
  10     FORMAT (a,'_',i4.4)
      ELSE
         WRITE (PRE_FILE_NAME,20) OUT_HISTORY%FNAME(1:LPOS)
  20     FORMAT (a)
      ENDIF
      IF (MULTIPLEOUT_HIS.eq.0) THEN
         WRITE (FILE_NAME,30) TRIM(PRE_FILE_NAME)
  30     FORMAT (a,'.nc')
      ELSE
# ifdef MPI_PARALL_GRID
         WRITE (FILE_NAME,40) TRIM(PRE_FILE_NAME),myrank
  40     FORMAT (a,'_',i4.4,'.nc')
# endif
      ENDIF
!
! Creating the initial history file description
!
      IF (IsInitDone.eqv..FALSE.) THEN ! At the beginning ...
        IsInitDone=.TRUE.
        recs_his2=0
!$OMP MASTER
        IF (WriteOutputProcess_his) THEN
! create nc file, vars, and do all time independant job
          iret = nf90_create(TRIM(FILE_NAME), NF90_CLOBBER, ncid)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
          nbTime=-1
          CALL WRITE_NETCDF_HEADERS_1(ncid, -1, MULTIPLEOUT_HIS, GRIDWRITE, IOBPD_HISTORY, np_write, ne_write)
          iret=nf90_inq_dimid(ncid, 'mnp', nnode_dims)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)
          iret=nf90_inq_dimid(ncid, 'ocean_time', ntime_dims)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)
          DO I=1,OUTVARS_COMPLETE
            IF (VAROUT_HISTORY%LVAR(I)) THEN
              CALL NAMEVARIABLE(I, eStr, eStrFullName, eStrUnit)
              iret=nf90_def_var(ncid,TRIM(eStr),NF90_OUTTYPE_HIS,(/ nnode_dims, ntime_dims /),var_id)
              CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
              iret=nf90_put_att(ncid,var_id,UNITS,TRIM(eStrUnit))
              CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
              iret=nf90_put_att(ncid,var_id,FULLNAME,TRIM(eStrFullName))
              CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
              IF (LSPHE) THEN
                iret=nf90_put_att(ncid,var_id,COORD,"lon lat")
              ELSE
                iret=nf90_put_att(ncid,var_id,COORD,"x y")
              END IF
              CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
            END IF
          END DO
          iret=nf90_close(ncid)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, iret)
        ENDIF
        CALL WRITE_NETCDF_HEADERS_2(FILE_NAME, MULTIPLEOUT_HIS, WriteOutputProcess_his, GRIDWRITE, np_write, ne_write)
!$OMP END MASTER
      END IF
!
! Computing the variables for output.
!
# ifdef MPI_PARALL_GRID
      IF (MULTIPLEOUT_HIS.eq.0) THEN
        ALLOCATE(OUTT_LOC(NP_GLOBAL,OUTVARS_COMPLETE), OUTT(NP_GLOBAL,OUTVARS_COMPLETE), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_output, allocate error 4')
        OUTT_LOC=0
        DO IP = 1, MNP
          ACLOC(:,:) = AC2(:,:,IP)
          CALL PAR_COMPLETE(IP, MSC, ACLOC, OUTPAR)
          OUTT_LOC(iplg(IP),:) = OUTPAR(:)
        END DO
        call mpi_reduce(OUTT_LOC,OUTT,NP_GLOBAL*OUTVARS_COMPLETE,rtype, MPI_SUM,0,comm,ierr)
        IF (myrank.eq.0) THEN
          DO IP=1,NP_GLOBAL
            OUTT(IP,:)=OUTT(IP,:)*nwild_gb(IP)
          enddo
        END IF
        DEALLOCATE(OUTT_LOC)
      ELSE
        ALLOCATE(OUTT(NP_RES,OUTVARS_COMPLETE), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_output, allocate error 5')
        DO IP = 1, NP_RES
          ACLOC(:,:) = AC2(:,:,IP)
          CALL PAR_COMPLETE(IP, MSC, ACLOC, OUTPAR)
          OUTT(IP,:) = OUTPAR(:)
        END DO
      ENDIF
# else
!$OMP PARALLEL DEFAULT(NONE) SHARED(AC2,OUTT,FRHIGH,MNP,MSC) PRIVATE(ACLOC,OUTPAR,IP)
!$OMP DO
      DO IP = 1, MNP
        ACLOC(:,:) = AC2(:,:,IP)
        CALL PAR_COMPLETE(IP, MSC, ACLOC, OUTPAR)
        OUTT(IP,:) = OUTPAR(:)
      END DO
!$OMP END DO
!$OMP END PARALLEL
# endif
      IF (LMONO_OUT) THEN
        OUTT(:,1) = OUTT(:,1) / SQRT(2.)
      ENDIF
!
! Writing down the variable in the history file.
!
      IF (IOBPD_HISTORY) THEN
        allocate(IOBPDoutput(MDC, np_write), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_output, allocate error 6')
        CALL GET_IOBPD_OUTPUT(IOBPDoutput, np_write)
      END IF
      recs_his2=recs_his2 + 1
      IF (WriteOutputProcess_his) THEN
!$OMP MASTER
        iret=nf90_open(TRIM(FILE_NAME), nf90_write, ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, iret)
        iret=nf90_inquire(ncid, unlimitedDimId = irec_dim)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, iret)
        iret=nf90_inquire_dimension(ncid, irec_dim,len = recs_his)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, iret)
        recs_his=recs_his+1
        IF (recs_his.ne.recs_his2) THEN
           CALL WWM_ABORT('There are more bugs to be solved');
        ENDIF
        CALL WRITE_NETCDF_TIME(ncid, recs_his, eTimeDay)
        IF (IOBPD_HISTORY) THEN
          iret=nf90_inq_varid(ncid, "IOBPD", var_id)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, iret)
          iret=nf90_put_var(ncid,var_id,IOBPDoutput,start = (/1, 1, recs_his/), count = (/ MDC, np_write, 1 /))
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, iret)
        END IF
        DO I=1,OUTVARS_COMPLETE
          IF (VAROUT_HISTORY%LVAR(I)) THEN
            CALL NAMEVARIABLE(I, eStr, eStrFullName, eStrUnit)
            iret=nf90_inq_varid(ncid, TRIM(eStr), var_id)
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 15, iret)
            IF (PRINTMMA) THEN
              CALL HISTORY_NC_PRINTMMA(eStr, OUTT,np_write,OUTVARS_COMPLETE,I)
            END IF
            IF (NF90_RUNTYPE == NF90_OUTTYPE_HIS) THEN
              iret=nf90_put_var(ncid,var_id,OUTT(:,I),start = (/1, recs_his/), count = (/ np_write, 1 /))
              CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 16, iret)
            ELSE
              iret=nf90_put_var(ncid,var_id,SNGL(OUTT(:,I)),start = (/1, recs_his/), count = (/ np_write, 1 /))
              CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 17, iret)
            ENDIF
          END IF
        END DO
        iret=nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 18, iret)
!$OMP END MASTER
      ENDIF
      IF (IOBPD_HISTORY) THEN
        deallocate(IOBPDoutput)
      END IF
# ifdef MPI_PARALL_GRID
      DEALLOCATE(OUTT)
# endif
      IF (OUT_HISTORY%IDEF.gt.0) THEN
        IF (recs_his2.eq.OUT_HISTORY%IDEF) THEN
          ifile=ifile+1
          IsInitDone = .FALSE.
        ENDIF
      ENDIF
      WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH OUTPUT_HISTORY_NC'
      FLUSH(STAT%FHNDL)
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
