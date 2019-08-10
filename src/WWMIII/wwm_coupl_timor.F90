#include "wwm_functions.h"
#if !defined ROMS_WWM_PGMCL_COUPLING && !defined MODEL_COUPLING_ATM_WAV && !defined MODEL_COUPLING_OCN_WAV
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_PIPES_TIMOR()
         USE DATAPOOL
         IMPLICIT NONE
!
! open pipe data files for coupling
!
           LSEWL = .TRUE.
           LSECU = .TRUE.
           WRITE(DBG%FHNDL,'("+TRACE...",A)') 'OPEN PIPE'
!          Pipes that are read by the wave model
           OPEN(1000,file='p_velx.dat'  ,form='unformatted', action='read')
           OPEN(1001,file='p_vely.dat'  ,form='unformatted', action='read')
           OPEN(1002,file='p_lev.dat'   ,form='unformatted', action='read')
           OPEN(1003,file='p_bot.dat'   ,form='unformatted', action='read')
           OPEN(1004,file='p_time.dat'  ,form='unformatted', action='write')
           OPEN(1005,file='p_dthyd.dat',form='unformatted', action='read')
!          Pipes that are written by the wave modell
           OPEN(101 ,file='p_stressx.dat' ,form='unformatted', action='write')
           OPEN(102 ,file='p_stressy.dat' ,form='unformatted', action='write')
           OPEN(103 ,file='p_waveh.dat'   ,form='unformatted', action='write')
           OPEN(104 ,file='p_wavet.dat'   ,form='unformatted', action='write')
           OPEN(105 ,file='p_waved.dat'   ,form='unformatted', action='write')
           OPEN(106 ,file='p_wavekm.dat'  ,form='unformatted', action='write')
           OPEN(107 ,file='p_wavetp.dat'  ,form='unformatted', action='write')
           OPEN(108 ,file='p_wavekp.dat'  ,form='unformatted', action='write')
           OPEN(109 ,file='p_orbit.dat'   ,form='unformatted', action='write')
           OPEN(110 ,file='p_stokesx.dat' ,form='unformatted', action='write')
           OPEN(111 ,file='p_stokesy.dat' ,form='unformatted', action='write')
           OPEN(112 ,file='p_windx.dat'   ,form='unformatted', action='write')
           OPEN(113 ,file='p_windy.dat'   ,form='unformatted', action='write')
!          Pipes writen as ergzus.bin for XFN
           OPEN(2003, file='pipe_wave.bin' ,form='unformatted' ,status = 'unknown')
           OPEN(2004, file='pipe_orbi.bin' ,form='unformatted' ,status = 'unknown')
           OPEN(2005, file='pipe_stok.bin' ,form='unformatted' ,status = 'unknown')
           WRITE(DBG%FHNDL,'("+TRACE...",A)') 'END OPEN PIPE'
!
! write coupling tyme step
!
           WRITE(1004) MAIN%DTCOUP
           FLUSH(1004)
!
! read coupling tyme step
!
           READ(1005)  MAIN%DTCUR

           IF(ABS(MAIN%DTCOUP-INT(MAIN%DTCOUP/MAIN%DTCUR)*MAIN%DTCUR).GT. .001) THEN
             write(DBG%FHNDL,*) 'TIME STEP OF THE hydraulic flow MODEL CANNOT BE DIVIDIED WITHOUT A REST'
             write(DBG%FHNDL,*)'dt Stroemung (s) =',MAIN%DTCUR, ',  dt Kopplung (s) = ',MAIN%DTCOUP
!             CALL WWM_ABORT('dt Stroemungsmodell muss gerades Vielfaches des Kopplungs-dt sein')
           END IF

           WRITE(DBG%FHNDL,'("+TRACE... DTCUR and DTCOUP",A)') MAIN%DTCUR, MAIN%DTCOUP

       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE TERMINATE_PIPES_TIMOR()
      USE DATAPOOL
      IMPLICIT NONE
      close(1000)
      close(1001)
      close(1002)
      close(1003)
      close(1004)
      close(1005)
      close(101)
      close(102)
      close(103)
      close(104)
      close(105)
      close(106)
      close(107)
      close(108)
      close(109)
      close(110)
      close(111)
      close(112)
      close(113)
      close(2003)
      close(2004)
      close(2005)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PIPE_TIMOR_IN(K)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: K
      IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN
        WATLEVOLD=WATLEV
        LCALC=.TRUE.
        WRITE(DBG%FHNDL,'("+TRACE...",A)') 'READING PIPE'
        WRITE(DBG%FHNDL,'("+TRACE...",A)') 'END READING PIPE'
        DEPDT = (WATLEV - WATLEVOLD) / MAIN%DTCOUP
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PIPE_TIMOR_OUT(K)
         IMPLICIT NONE

         INTEGER, INTENT(IN)  :: K

      END SUBROUTINE
#endif