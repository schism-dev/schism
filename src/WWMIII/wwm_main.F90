#include "wwm_functions.h"
#include "wwm_functions.h"
!:)
! __      __  __      __  _____  .___.___.___ 
!/  \    /  \/  \    /  \/     \ |   |   |   |
!\   \/\/   /\   \/\/   /  \ /  \|   |   |   |
! \        /  \        /    Y    \   |   |   |
!  \__/\  /    \__/\  /\____|__  /___|___|___|
!       \/          \/         \/             
!
! WWM-III (Wind Wave Model) source code 
! 
! The 1st version of the WWM code was written by Jian-Ming Liau in his thesis supervised by Tai-Wen Hsu (Liau et al. 2002). 
! The source code served as the basis for my thesis that was as well supervised by Tai-Wen Hsu and Ulrich Zanke. In my thesis 
! new numerics and source terms have beend developed (Roland, 2008) and resulted in the WWM-II version of the code. Following this
! the code has served from than as a basis for a 10 year development. In this time the source code was significantly rewritten and 
! enhanced with various capabilities. The numerics have been completely revised (Roland, 2008) for WWM-II.
!
! In WWM-III the biggest contributions was the work with Joseph Zhang, involving fully coupling to SCHISM via fine data structures and 
! introducing domain decomposition techniques with the SCHISM coupling and later together with Thomas Huxhorn via PDLIB for WWM-III only. 
! Moreover, we have obtained severall source term packages by the original authors, which are widely used in spectral wave models such as 
! SWAN, WAM, or WW3. 
!
! the source term package of Ardhuin et al. 2009, 2010 (ST4, courtesy Fabrice Ardhuin) from ECWAM (courtesy Jean-Bidlot, Peter Janssen)i as 
! well as the ST6 observation based package by Stefan Zieger, Zieger et al. (2015) based on the work of Babanin et al. (XXXX, XXXX, XXXX, XXXX, XXXX)
! was implemented in the WWM-III. 
!
! The I/O was completely rewritten in NETCDF and various common wind fields can be read such as CFRS, ECMWF, NCEP or others.  
! Parallelization is done using the PDLIB decomposition library developed by Roland & Partner and based on domain decmoposition using ParMetis. 
! 
! Developers: The WWM Group
!
! Lead:  
! -----
! Jian-Ming Liau (XXXX, Taiwan, WWM-I)
! Aron Roland (Roland & Partner, Darmstadt, WWM-II, WWM-III),
! Jian-Ming Liau (XXXX, Taiwan, WWM-I)
! Mathieu Dutour Sikiric (IRB, Zagreb, WWM-III),
! Yinglong Joseph Zhang (VIMS, USA, WWM-III),
!
! Contributors:
! -------------
! Fabrice Ardhuin (IFREMER, Brest, France),
! Alexander Babanin (XXXX, Australia) 
! Xavier Bertin (UNR, La Rochelle, France), 
! Jean Bidlot (ECMWF, Reading, U.K.)
! Guillaume Dodet (IFREMER, Brest, France),
! Christian Ferrarin (ISMAR-CNR, Venice, Italy),
! Andrea Fortunato (LNEC, Lissabon, Portugal),
! Thomas Huxhorn (BGS IT&E, Darmstadt, Germany),
! Ivica Janekovic, (UWA, Perth, Australia) 
! Kai Li (XXXX, PR China), 
! Kevin Maarten (UNR, La Rochele, France), 
! Peter Janssen (ECWMF, Reading, U.K.) 
! Stefan Zieger (XXXX, Australia) 
!
! Initiators:
! -----------
! Tai-Wen Hsu (NTOU, Taiwan, WWM-I) 
! Ulrich Zanke (Zanke & Partner, Garbsen, Germany, WWM-II)
!       
! Copyright: 2008 - 2019
! All Rights Reserved                                     
!
! http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef SCHISM
 !!!     SUBROUTINE WWM_II(IT_SCHISM,icou_elfe_wwm,DT_SELFE0,NSTEP_WWM0)
      SUBROUTINE WWM_II(IT_SCHISM,icou_elfe_wwm,DT_SCHISM0,NSTEP_WWM0,RADFLAG2)

         USE DATAPOOL
         use  schism_msgp !, only : myrank,parallel_abort,itype,comm,ierr
         use schism_glbl, only : iplg,wwave_force

         IMPLICIT NONE

         INTEGER, INTENT(IN)           :: NSTEP_WWM0, icou_elfe_wwm
         REAL(rkind), INTENT(IN)       :: DT_SCHISM0
         CHARACTER(LEN=3), INTENT(OUT) :: RADFLAG2

         REAL(rkind), SAVE  :: SIMUTIME
         REAL(rkind)        :: T1, T2
         REAL(rkind)        :: TIME1, TIME2, TIME3, TIME4, TIME5, TIME6, TIME7

         INTEGER     :: I, IP, IT_SCHISM, K, IFILE, IT
         REAL(rkind) :: DT_PROVIDED
         REAL(rkind) :: OUTPAR(OUTVARS), OUTWINDPAR(WINDVARS), ACLOC(MSC,MDC)
         character(LEN=15) :: CALLFROM

         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING WWM_II'
         CALL FLUSH(STAT%FHNDL)

#ifdef TIMINGS
         TIME1 = mpi_wtime()
#endif 

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' STARTING WWM FROM SCHISM ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) call wwm_abort('NAN IN MAIN 1')
         ENDIF

         if(WINDVARS/=size(WIND_INTPAR,2)) call wwm_abort('Dimension mismatch: OUTWINDPAR and out_wwm_windpar')
         if(OUTVARS/=size(OUTT_INTPAR,2))  call wwm_abort('Dimension mismatch: OUTPAR and out_wwm')

         NSTEPWWM = NSTEP_WWM0

         DT_SCHISM      = DT_SCHISM0

#ifdef TIMINGS
         T1 = MyREAL(IT_SCHISM-NSTEPWWM)*DT_SELFE0 ! Beginn time step ...
         T2 = MyREAL(IT_SCHISM)*DT_SELFE0          ! End of time time step ...
#endif 

         DT_PROVIDED=NSTEPWWM*DT_SCHISM

         IF (abs(MAIN%DELT - DT_PROVIDED).gt.THR) THEN
           WRITE(DBG%FHNDL,*) 'MAIN%DELT=', MAIN%DELT, ' in wwminput.nml'
           WRITE(DBG%FHNDL,*) 'But nstep_wwm*dt=', DT_PROVIDED
           WRITE(DBG%FHNDL,*) 'nstep_wwm=', NSTEPWWM
           WRITE(DBG%FHNDL,*) '       dt=', DT_SCHISM
           CALL WWM_ABORT('Correct coupled model time-steppings')
         ENDIF

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' FIRST SUM IN MAIN ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN MAIN 2')
         ENDIF

         WRITE(STAT%FHNDL,'("+TRACE...",A)') ' ---- ALL CHECKS DONE'
         CALL FLUSH(STAT%FHNDL)

         SIMUTIME = SIMUTIME + MAIN%DELT

         IF (icou_elfe_wwm == 1) THEN ! Full coupling 
           WLDEP       = DEP8
           WATLEV      = ETA2
           WATLEVOLD   = ETA1
           DEP         = MAX(ZERO,WLDEP + WATLEV)
           CURTXY(:,1) = UU2(NVRT,:)
           CURTXY(:,2) = VV2(NVRT,:)
           IF (.NOT. LWINDFROMWWM) THEN
             WINDXY(:,1) = WINDX0
             WINDXY(:,2) = WINDY0
           END IF
           LSECU       = .TRUE.
           LSEWL       = .TRUE.
           LCALC       = .TRUE.
         ELSE IF (icou_elfe_wwm == 0) THEN ! No interaction at all 
           WLDEP       = DEP8
           WATLEV      = ZERO 
           WATLEVOLD   = ZERO
           DEP         = WLDEP
           CURTXY(:,1) = ZERO !REAL(rkind)(UU2(NVRT,:))
           CURTXY(:,2) = ZERO !REAL(rkind)(VV2(NVRT,:))
           IF (.NOT. LWINDFROMWWM) THEN
             WINDXY(:,1) = WINDX0
             WINDXY(:,2) = WINDY0
           END IF
           LSECU       = .FALSE.
           LSEWL       = .FALSE.
           LCALC       = .TRUE. 
         ELSE IF (icou_elfe_wwm == 2) THEN ! Currents and water levels in wwm but no radiation stress in SCHISM
           WLDEP       = DEP8
           WATLEV      = ETA2
           WATLEVOLD   = ETA1
           DEP         = MAX(ZERO, WLDEP + WATLEV)
           CURTXY(:,1) = UU2(NVRT,:)
           CURTXY(:,2) = VV2(NVRT,:)
           IF (.NOT. LWINDFROMWWM) THEN
             WINDXY(:,1) = WINDX0
             WINDXY(:,2) = WINDY0
           END IF
           LSECU       = .TRUE.
           LSEWL       = .TRUE.
           LCALC       = .TRUE.
         ELSE IF (icou_elfe_wwm == 3) THEN ! No current and no water levels in wwm but radiation stress in SCHISM
           WLDEP       = DEP8
           WATLEV      = ZERO
           WATLEVOLD   = ZERO
           DEP         = WLDEP
           CURTXY(:,1) = ZERO !REAL(rkind)(UU2(NVRT,:))
           CURTXY(:,2) = ZERO !REAL(rkind)(VV2(NVRT,:))
           IF (.NOT. LWINDFROMWWM) THEN
             WINDXY(:,1) = WINDX0
             WINDXY(:,2) = WINDY0
           END IF
           LSECU       = .FALSE.
           LSEWL       = .FALSE.
           LCALC       = .TRUE.
         ELSE IF (icou_elfe_wwm == 4 .or. icou_elfe_wwm == 5) THEN ! No current but water levels in wwm and radiation stresss in SCHISM
           WLDEP       = DEP8
           WATLEV      = ETA2
           WATLEVOLD   = ETA1
           DEP         = WLDEP
           CURTXY(:,1) = 0.!UU2(NVRT,:) 
           CURTXY(:,2) = 0.!UU2(NVRT,:) 
           IF (.NOT. LWINDFROMWWM) THEN
             WINDXY(:,1) = WINDX0
             WINDXY(:,2) = WINDY0
           END IF
           LSECU       = .FALSE.
           LSEWL       = .TRUE.
           LCALC       = .TRUE.
         ELSE IF (icou_elfe_wwm == 6) THEN ! Currents but no water levels in wwm and radiation stress in SCHISM  
           WLDEP       = DEP
           WATLEV      = ZERO 
           WATLEVOLD   = ZERO 
           DEP         = WLDEP
           CURTXY(:,1) = UU2(NVRT,:) 
           CURTXY(:,2) = VV2(NVRT,:) 
           IF (.NOT. LWINDFROMWWM) THEN
             WINDXY(:,1) = WINDX0
             WINDXY(:,2) = WINDY0
           END IF
           LSECU       = .TRUE.
           LSEWL       = .FALSE.
           LCALC       = .TRUE.
         ELSE IF (icou_elfe_wwm == 7) THEN ! Currents but no water levels in wwm and no radiation stress in SCHISM  
           WLDEP       = DEP
           WATLEV      = ZERO
           WATLEVOLD   = ZERO
           DEP         = WLDEP
           CURTXY(:,1) = UU2(NVRT,:)
           CURTXY(:,2) = UU2(NVRT,:)
           IF (.NOT. LWINDFROMWWM) THEN
             WINDXY(:,1) = WINDX0
             WINDXY(:,2) = WINDY0
           END IF
           LSECU       = .TRUE.
           LSEWL       = .FALSE.
           LCALC       = .TRUE.
         END IF

         DEPDT = (WATLEV - WATLEVOLD) / DT_SCHISM0

         IF (LNANINFCHK) THEN
           CALL SCHISM_NANCHECK_INPUT_A
         END IF

         IFILE = 1
         IT    = 1
         IF (LBCSE) THEN
           CALL SET_WAVE_BOUNDARY_CONDITION
         END IF

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER SETTING BOUNDARY CONDITION IN MAIN ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN MAIN 3')
         ENDIF

         IF (LFIRSTSTEP) THEN
           IF (INITSTYLE == 1) CALL INITIAL_CONDITION!We need to call for the case of wind dependent intial guess this call since before we have no wind from SCHISM
           LFIRSTSTEP = .FALSE.
           LCALC      = .TRUE.
         END IF

#ifdef TIMINGS
         TIME2 = mpi_wtime() 
#endif

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' BEFORE COMPUTE ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN MAIN 4')
         ENDIF

         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE'
         CALL FLUSH(STAT%FHNDL)

         CALLFROM='SCHISM'
         IF (LQSTEA) THEN
           CALL QUASI_STEADY(KKK)
         ELSE
           CALL UN_STEADY(KKK,CALLFROM)
         END IF

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER COMPUTE ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN MAIN 5') 
         ENDIF

#ifdef TIMINGS
         TIME3 = mpi_wtime()
#endif 

         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'FINISHED COMPUTE nth call to WWM', SIMUTIME
         CALL FLUSH(STAT%FHNDL)

         DO IP = 1, MNP
           ACLOC = AC2(:,:,IP)
           IF (DEP(IP) .GT. DMIN) THEN
             CALL INTPAR(IP, MSC, ACLOC, OUTPAR)
             OUTT_INTPAR(IP,:) = OUTPAR
             CALL WINDPAR(IP,OUTWINDPAR)
             WIND_INTPAR(IP,:) = OUTWINDPAR
             IF (LMONO_OUT) THEN
               OUTT_INTPAR(IP,1) = OUTT_INTPAR(IP,1) / SQRT(TWO)
             END IF
           ELSE
             OUTT_INTPAR(IP,:) = ZERO
             WIND_INTPAR(IP,:) = ZERO
           END IF
         END DO

         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'FINISHED FILLING RESULTS', SIMUTIME
         CALL FLUSH(STAT%FHNDL)

#ifdef TIMINGS
         TIME4 = mpi_wtime()
#endif

!
! Compute wave forces ...
!
         ! Initialization
         RADFLAG2 = RADFLAG !for output into SCHISM
!         WWAVE_FORCE = ZERO
         IF (icou_elfe_wwm == 0 .OR. icou_elfe_wwm == 2 .OR. icou_elfe_wwm == 5 .OR. icou_elfe_wwm == 7) THEN
           STOKES_VEL = ZERO
           JPRESS = ZERO
           SBR = ZERO
           SBF = ZERO
           !STOKES_X=ZERO
           !STOKES_Y=ZERO
         ELSE 
           IF (RADFLAG == 'VOR') THEN                  ! Vortex force formalism as described in Bennis (2011)
             CALL STOKES_STRESS_INTEGRAL_SCHISM        ! Compute Stokes drift velocities and pressure terms 
             CALL COMPUTE_CONSERVATIVE_VF_TERMS_SCHISM ! Conservative terms (relative to Stokes drift advection, Coriolis and pressure head: Eq. 17, 19 and 20 from Bennis 2011)
             CALL COMPUTE_BREAKING_VF_TERMS_SCHISM     ! Sink of momentum due to wave breaking and update wwave_force
           ELSE ! Radiation stress formalism (Longuet-Higgins and Stewart, 1962 and 1964) as described in Battjes (1974)
             CALL RADIATION_STRESS_SCHISM
           ENDIF
		   
           !...  Fixing the wave forces at the shoreline (dry/wet boundary)
           !     If shorewafo == 1, we impose the wave forces to be equal to the  barotropic gradient
           !     The wwave_force values at the shoreline computed above with hgrad_nodes are hence overwritten
           IF(SHOREWAFO == 1) CALL SHORELINE_WAVE_FORCES
		   
           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'FINISHED FILLING VORTEX', SIMUTIME
           CALL FLUSH(STAT%FHNDL)
         END IF 

#ifdef TIMINGS
         TIME5 = mpi_wtime()
#endif
         IF (LNANINFCHK) THEN
           CALL SCHISM_NANCHECK_INPUT_B
         END IF

         KKK = KKK + 1

#ifdef TIMINGS
         TIME6 = mpi_wtime()
#endif

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' END OF MAIN ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT ('NAN IN MAIN 5')
         ENDIF

         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'END OF COMPUTATIONS NOW RETURN TO SCHISM', SIMUTIME
         CALL FLUSH(STAT%FHNDL)

#ifdef TIMINGS
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----TOTAL TIMINGS-----'
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'PREPARATION        ', TIME2-TIME1
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'INTEGRATION        ', TIME3-TIME2
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'OUTPUT TO SCHISM   ', TIME4-TIME3
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'RADIATION STRESSES ', TIME5-TIME4
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'NAN CHECK          ', TIME6-TIME5
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'TOTAL TIME         ', TIME6-TIME1
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '------END-TIMINGS-  ---'
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'FINISHED WITH WWM', SIMUTIME
         CALL FLUSH(STAT%FHNDL)
#endif
 
      END SUBROUTINE WWM_II
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SCHISM_NANCHECK_INPUT_A
      USE DATAPOOL
      implicit none
      integer IP
      DO IP = 1, MNP
        IF (WINDXY(IP,1) .NE. WINDXY(IP,1)) THEN
          WRITE(DBG%FHNDL,*) 'NaN in WINDX', IP, WINDXY(IP,1) 
          CALL FLUSH(DBG%FHNDL)
        END IF
        IF (WINDXY(IP,2) .NE. WINDXY(IP,2)) THEN
          WRITE(DBG%FHNDL,*) 'NaN in WINDY', IP, WINDXY(IP,2) 
          CALL FLUSH(DBG%FHNDL)
        END IF
        IF (WATLEV(IP) .NE. WATLEV(IP)) THEN
          WRITE(DBG%FHNDL,*) 'NaN in WATLEV', IP, WATLEV(IP) 
          CALL FLUSH(DBG%FHNDL)
        END IF
        IF (WATLEVOLD(IP) .NE. WATLEVOLD(IP)) THEN
          WRITE(DBG%FHNDL,*) 'NaN in WATLEV', IP, WATLEV(IP)
          CALL FLUSH(DBG%FHNDL)
        END IF
        IF (CURTXY(IP,1) .NE. CURTXY(IP,1)) THEN
          WRITE(DBG%FHNDL,*) 'NaN in CURTX', IP, CURTXY(IP,1)
          CALL FLUSH(DBG%FHNDL)
        END IF
        IF (CURTXY(IP,2) .NE. CURTXY(IP,2)) THEN
          WRITE(DBG%FHNDL,*) 'NaN in CURTY', IP, CURTXY(IP,2)
          CALL FLUSH(DBG%FHNDL)
        END IF
      END DO
      END SUBROUTINE SCHISM_NANCHECK_INPUT_A
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SCHISM_NANCHECK_INPUT_B
      USE DATAPOOL
      implicit none
      integer IP, I
      DO IP = 1, MNP
        IF (SUM(OUTT_INTPAR(IP,:)) .NE. SUM(OUTT_INTPAR(IP,:))) THEN
          DO I = 1, SIZE(OUTT_INTPAR(IP,:))
            WRITE(DBG%FHNDL,*) 'NaN in OUTT_INTPAR', IP, I, OUTT_INTPAR(IP,I)
            CALL FLUSH(DBG%FHNDL)
          END DO
        END IF
        IF (SUM(WIND_INTPAR(IP,:)) .NE. SUM(WIND_INTPAR(IP,:))) THEN
          DO I = 1, SIZE(WIND_INTPAR(IP,:))
            WRITE(DBG%FHNDL,*) 'NaN in WIND_INTPAR', IP, I, WIND_INTPAR(IP,I)
            CALL FLUSH(DBG%FHNDL)
          END DO
        END IF
!Error: force defined at side centers
        IF (SUM(WWAVE_FORCE(:,IP,:)) .NE. SUM(WWAVE_FORCE(:,IP,:))) THEN
          DO I = 1, SIZE(WWAVE_FORCE(1,:,IP)) ! loop over layers ...
            WRITE(DBG%FHNDL,*) 'NaN in WWAVE_FORCE', IP, I, WWAVE_FORCE(1,I,IP), WWAVE_FORCE(2,I,IP)
            CALL FLUSH(DBG%FHNDL)
          END DO
        END IF 
      END DO
      END SUBROUTINE SCHISM_NANCHECK_INPUT_B
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE UN_STEADY(K,CALLFROM)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: K

         REAL(rkind)    :: CONV1, CONV2, CONV3, CONV4, CONV5
         REAL(rkind)    :: TIME1, TIME2, TIME3, TIME4, TIME5, TIME6

         CHARACTER(LEN=15)   :: CTIME,CALLFROM

#ifdef TIMINGS
      CALL WAV_MY_WTIME(TIME1)
#endif

      CALL IO_1(K)

#ifdef TIMINGS
      CALL WAV_MY_WTIME(TIME2)
#endif

!      CALL Print_SumAC2("Before the advection")
      IF (ICOMP .EQ. 0) THEN
        CALL COMPUTE_SIMPLE_EXPLICIT
      ELSE IF (ICOMP .EQ. 1) THEN 
        CALL COMPUTE_SEMI_IMPLICIT
      ELSE IF (ICOMP .EQ. 2) THEN 
        CALL COMPUTE_SEMI_IMPLICIT
      ELSE IF (ICOMP .EQ. 3) THEN 
        CALL COMPUTE_IMPLICIT
      END IF
!      CALL Print_SumAC2("After the advection")

#ifdef TIMINGS
      CALL WAV_MY_WTIME(TIME3)
#endif

#ifdef WWM_SETUP
      IF (LZETA_SETUP) THEN
        CALL WAVE_SETUP_COMPUTATION
      END IF
#endif

#ifdef TIMINGS
      CALL WAV_MY_WTIME(TIME4)
#endif

      MAIN%TMJD = MAIN%BMJD + MyREAL(K)*MAIN%DELT*SEC2DAY
      RTIME = MAIN%TMJD - MAIN%BMJD

#ifndef SCHISM
#if defined WWM_MPI
      IF (myrank.eq.0) WRITE(*,101)  K, MAIN%ISTP, RTIME
#else
      WRITE(12,101)  K, MAIN%ISTP, RTIME
#endif 
#endif
      CALL IO_2(K)
!      CALL Print_SumAC2("After IO_2")

#ifdef TIMINGS
      CALL WAV_MY_WTIME(TIME5)
#endif

      IF (LCONV) THEN
        CALL CHECK_STEADY(RTIME,CONV1,CONV2,CONV3,CONV4,CONV5)
      END IF

#ifdef TIMINGS
      CALL WAV_MY_WTIME(TIME6)
#endif

      IF (.NOT. LDIFR) LCALC = .FALSE.

      CALL MJD2CT(MAIN%TMJD, CTIME)

#ifdef TIMINGS
# ifdef MPI_PARALL_GRID
      IF (myrank == 0) THEN
# endif
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6,A20)') '-----SIMULATION TIME-----        ', MAIN%TMJD, CTIME
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----TOTAL RUN TIMES-----        '
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'PREPROCESSING                    ', TIME2-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'INTEGRATION                      ', TIME3-TIME2
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'WAVE SETUP                       ', TIME4-TIME3
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'POSTPROCESSING                   ', TIME5-TIME4
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CHECK STEADY                     ', TIME6-TIME5
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'TOTAL TIME                       ', TIME6-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-------------TIMINGS-------------'
        FLUSH(STAT%FHNDL)
# ifdef MPI_PARALL_GRID
      ENDIF
# endif
#endif

      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'LEAVING UN_STEADY'

101      FORMAT ('+STEP = ',I10,'/',I10,' ( TIME = ',F15.4,' DAYS)')
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE QUASI_STEADY(K)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: K
      INTEGER :: IT
      REAL(rkind)    :: ITERTIME
      REAL(rkind)    :: CONV1, CONV2, CONV3, CONV4, CONV5

      CALL IO_1(K)

      IF (LCFL) CALL CFLSPEC()

      IF (LCHKCONV) IP_IS_STEADY = 0 ! Reset local convergence indicators ...
      IF (LCHKCONV) IE_IS_STEADY = 0

#ifdef MPI_PARALL_GRID
!         NQSITER = NSTEPWWM ! this is not very flexible!
#endif

      DO IT = 1, NQSITER

        DT_ITER = MAIN%DELT/MyREAL(NQSITER)

        IF (ICOMP .EQ. 0) THEN
          CALL COMPUTE_SIMPLE_EXPLICIT
        ELSE IF (ICOMP .EQ. 1) THEN
          CALL COMPUTE_SEMI_IMPLICIT
        ELSE IF (ICOMP .EQ. 2) THEN
          CALL COMPUTE_SEMI_IMPLICIT
        ELSE IF (ICOMP .EQ. 3) THEN
          CALL COMPUTE_IMPLICIT
        END IF

        ITERTIME = RTIME*DAY2SEC+IT*DT_ITER

        IF (LCHKCONV) THEN
          CALL CHECK_STEADY(ITERTIME,CONV1,CONV2,CONV3,CONV4,CONV5)
!             DO IP = 1, MNP
!               IF (IP_IS_STEADY(IP) .GE. 1) AC2(:,:,IP) = AC1(IP,:,:)
!             ENDDO
          IF ( (CONV1 .GT. 100._rkind*QSCONV1 .AND.                       &
     &             CONV2 .GT. 100._rkind*QSCONV2 .AND.                    &
     &             CONV3 .GT. 100._rkind*QSCONV3 .AND.                    &
     &             CONV4 .GT. 100._rkind*QSCONV4 .AND.                    &
     &             CONV5 .GT. 100._rkind*QSCONV5 .AND.                    &
     &             K .NE. 1) .OR.                                         &
     &             IT .EQ. NQSITER ) THEN
#ifndef SCHISM
            WRITE(QSTEA%FHNDL,'(3I10,5F15.8)') K, IT, NQSITER, CONV1, CONV2, CONV3, CONV4, CONV5
#else
            if (myrank == 0) WRITE(QSTEA%FHNDL,'(3I10,5F15.8)') K, IT, NQSITER, CONV1, CONV2, CONV3, CONV4, CONV5
#endif
            FLUSH(QSTEA%FHNDL)
            EXIT 
          END IF
        END IF
        IF (LOUTITER) CALL WWM_OUTPUT(ITERTIME,.FALSE.)
      END DO

#ifdef MPI_PARALL_GRID
      MAIN%TMJD = MAIN%BMJD + MyREAL(K)*MAIN%DELT*SEC2DAY
      RTIME = MAIN%TMJD - MAIN%BMJD
      IF (myrank == 0) WRITE(STAT%FHNDL,101)  K, MAIN%ISTP, RTIME*DAY2SEC
#else
      MAIN%TMJD = MAIN%BMJD + MyREAL(K)*MAIN%DELT*SEC2DAY
      RTIME = MAIN%TMJD - MAIN%BMJD
      WRITE(STAT%FHNDL,101)  K, MAIN%ISTP, RTIME*DAY2SEC
#endif
      FLUSH(STAT%FHNDL)
      CALL IO_2(K)
      IF (.NOT. LDIFR) LCALC = .FALSE.

101   FORMAT ('+STEP = ',I5,'/',I5,' ( TIME = ',F15.4,'HR )')
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE IO_1(K)
#if defined ROMS_WWM_PGMCL_COUPLING || defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
      USE WWMaOCN_PGMCL
#endif
#if defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
      USE pgmcl_lib_WWM, only : WAV_all_import_export
#endif
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: K
      INTEGER             :: IT, IFILE
      IF (LWINDFROMWWM) THEN
        CALL UPDATE_WIND(K)
      END IF
#ifndef SCHISM
      IF (.NOT. LCPL) THEN
        IF (LSECU) THEN
          CALL UPDATE_CURRENT(K)
        END IF
        IF (LSEWL) THEN
          CALL UPDATE_WATLEV(K)
        END IF
      END IF
#endif
#ifndef SCHISM
      IF (LBCSE) THEN
        CALL SET_WAVE_BOUNDARY_CONDITION
      END IF
#endif
!
!      *** coupling via pipe *** read pipe
!
#if !defined SCHISM && !defined ROMS_WWM_PGMCL_COUPLING && !defined MODEL_COUPLING_ATM_WAV && !defined MODEL_COUPLING_OCN_WAV
      IF (LCPL .AND. LTIMOR) THEN
        CALL PIPE_TIMOR_IN(K)
# ifdef SHYFEM_COUPLING
      ELSE IF (LCPL .AND. LSHYFEM) THEN
        CALL PIPE_SHYFEM_IN(K)
# endif
      ELSE IF (LCPL .AND. LROMS) THEN
        CALL PIPE_ROMS_IN(K,IFILE,IT)
      END IF
#endif
#ifdef ROMS_WWM_PGMCL_COUPLING
      IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN
        CALL WAV_ocnAwav_import(K,IFILE,IT)
      END IF
      IF (K == 1) CALL INITIAL_CONDITION
#endif
#if defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
      CALL WAV_all_import_export(K, IFILE, IT)
#endif
!
!      *** recalculate water level and current related values 
!
      IF (LSEWL .OR. LSECU .OR. LCPL) THEN ! LCPL makes sure that when the model is coupled it gets into this part for 100%
        DEP  = MAX(ZERO,WLDEP + WATLEV) ! d = a + h  if -h .gt. a set d to zero
        CALL SETSHALLOW
        CALL GRADDEP
        IF (MESTR == 6) CALL GRAD_CG_K
        CALL WAVE_K_C_CG
        CALL GRADCURT
        CALL SET_IOBPD
        CALL SET_IOBPD_BY_DEP
        IF (LCFL) THEN
          CFLCXY = ZERO
          CALL CFLSPEC
        ENDIF
        IF (LMAXETOT .AND. MESBR == 0) CALL SET_HMAX
      END IF
!
!
!
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE IO_2(K)
#if defined ROMS_WWM_PGMCL_COUPLING || defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
      USE WWMaOCN_PGMCL
#endif
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: K
      CALL GENERAL_OUTPUT
#ifndef SCHISM
# if !defined ROMS_WWM_PGMCL_COUPLING && !defined MODEL_COUPLING_ATM_WAV && !defined MODEL_COUPLING_OCN_WAV
      IF (LCPL .AND. LTIMOR) THEN
        CALL PIPE_TIMOR_OUT(K)
#  ifdef SHYFEM_COUPLING
      ELSE IF (LCPL .AND. LSHYFEM) THEN
        CALL PIPE_SHYFEM_OUT(K)
#  endif
      ELSE IF (LCPL .AND. LROMS) THEN
        CALL PIPE_ROMS_OUT(K)
      END IF
# endif
# ifdef ROMS_WWM_PGMCL_COUPLING
      IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN
        CALL WAV_ocnAwav_export(K)
      END IF
# endif
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#if !defined SCHISM && !defined PDLIB && defined MPI_PARALL_GRID
      SUBROUTINE SIMPLE_PRE_READ
      USE DATAPOOL
      USE schism_glbl, only : msc2, mdc2, ics
      IMPLICIT NONE
      CHARACTER(LEN=20) :: BEGTC, UNITC, ENDTC
      REAL(rkind) DELTC
         NAMELIST /PROC/ PROCNAME, DIMMODE, LSTEA, LQSTEA, LSPHE,       &
     &      LNAUTIN, LNAUTOUT, LMONO_OUT, LMONO_IN,                     &
     &      BEGTC, DELTC, UNITC, ENDTC, DMIN 
         NAMELIST /GRID/ LCIRD, LSTAG, MINDIR, MAXDIR, MDC, FRLOW,      &
     &      FRHIGH, MSC, FILEGRID, IGRIDTYPE, LSLOP, SLMAX, LVAR1D,     &
     &      LOPTSIG, CART2LATLON, LATLON2CART 
      INTEGER FHNDL
      !
      FHNDL=12
      CALL TEST_FILE_EXIST_DIE("Missing input file : ", TRIM(INP%FNAME))
      OPEN(FHNDL, FILE = TRIM(INP%FNAME))
      READ(FHNDL, NML = PROC)
      READ(FHNDL, NML = GRID)
      CLOSE(FHNDL)
      IF (LSPHE) THEN
        ics=2
      ELSE
        ics=1
      ENDIF
      IF (CART2LATLON) THEN
        ics=1
      END IF
      IF (LATLON2CART) THEN
        ics=2
      END IF
      IF (CART2LATLON .and. LATLON2CART) THEN
        CALL WWM_ABORT('You cannot have both CART2LATLON and CART2LONLAT')
      END IF
      msc2=MSC
      mdc2=MDC
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
#if !defined SCHISM
# if defined ROMS_WWM_PGMCL_COUPLING || defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
      SUBROUTINE WWMIII_MPI(MyCOMM)
# else
      PROGRAM WWMIII_MPI
# endif

# ifdef ROMS_WWM_PGMCL_COUPLING
      USE mod_coupler, only : WAV_COMM_WORLD, MyRankGlobal
# endif
# if defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
      USE coupling_var, only : WAV_COMM_WORLD, MyRankGlobal
# endif

      USE DATAPOOL, only: MAIN, SEBO,                                  &
     &      NDT_BND_FILE, IWBNDLC, AC2, WBAC, STAT, RTIME,             &
     &      bnd_time_all_files, LSPHE, WLDEP, DEP, SMALL, KKK,         &
     &      WATLEV, LBCSE, LBCWA, LBCSP, IWBMNP, IWBNDLC, WBAC,        &
     &      WBACOLD, WBACNEW, DSPEC, LBINTER, LFIRSTSTEP, LQSTEA,      &
     &      LINHOM, IBOUNDFORMAT, DAY2SEC, SEC2DAY,                    &
     &      NUM_NETCDF_FILES_BND, LSECU, RKIND, MDC, MSC, MNP

# ifdef MPI_PARALL_GRID
      USE datapool, only: rkind, comm, myrank, ierr, nproc,            &
     &      parallel_finalize
# endif

      implicit none

# if defined MPI_PARALL_GRID || defined ROMS_WWM_PGMCL_COUPLING || defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
      include 'mpif.h'
# endif
# if defined ROMS_WWM_PGMCL_COUPLING || defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
      integer, intent(in) :: MyCOMM
# endif
# ifdef TIMINGS 
      REAL(rkind)        :: TIME1, TIME2
# endif
      integer :: i,j,k, IP
      character(len=15) CALLFROM
# if defined DEBUG && (defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV)
      write(740+MyRankGlobal,*)  'WWMIII_MPI, before mpi_init'
      FLUSH(740+MyRankGlobal)
# endif

# if defined WWM_MPI && !defined ROMS_WWM_PGMCL_COUPLING && !defined MODEL_COUPLING_ATM_WAV && !defined MODEL_COUPLING_OCN_WAV
      call mpi_init(ierr)
      if(ierr/=MPI_SUCCESS) call wwm_abort('Error at mpi_init')
# endif
# if defined DEBUG && (defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV)
      write(740+MyRankGlobal,*)  'WWMIII_MPI, after mpi_init'
      FLUSH(740+MyRankGlobal)
# endif

# ifdef TIMINGS
      CALL WAV_MY_WTIME(TIME1)
# endif
# if defined DEBUG && (defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV)
      write(740+MyRankGlobal,*)  'WWMIII_MPI, after WAV_MY_WTIME'
      FLUSH(740+MyRankGlobal)
# endif
# if defined ROMS_WWM_PGMCL_COUPLING || defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
      comm=MyCOMM
      WAV_COMM_WORLD=MyCOMM
# else
#  ifdef WWM_MPI
      call mpi_comm_dup(MPI_COMM_WORLD,comm,ierr)
      if(ierr/=MPI_SUCCESS) call wwm_abort('Error at mpi_comm_dup')
#  endif
# endif
# if defined DEBUG && (defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV)
      write(740+MyRankGlobal,*)  'WWMIII_MPI, after mpi_comm_dup and WAV_COMM_WORLD'
      FLUSH(740+MyRankGlobal)
# endif
# ifdef MPI_PARALL_GRID
      call mpi_comm_size(comm,nproc,ierr)
      if(ierr/=MPI_SUCCESS) call wwm_abort('Error at mpi_comm_size')
      call mpi_comm_rank(comm,myrank,ierr)
      if(ierr/=MPI_SUCCESS) call wwm_abort('Error at mpi_comm_rank')
#  ifndef PDLIB 
      CALL SIMPLE_PRE_READ
#  endif
      CALLFROM='WWM_MPI'
# else
      CALLFROM='WWM'
# endif
# if defined DEBUG && (defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV)
      write(740+MyRankGlobal,*)  'WWMIII_MPI, after mpi_comm_size/rank'
      FLUSH(740+MyRankGlobal)
# endif
      CALL INITIALIZE_WWM
# if defined DEBUG && (defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV)
      write(740+MyRankGlobal,*)  'WWMIII_MPI, after INITIALIZE_WWM'
      FLUSH(740+MyRankGlobal)
# endif

!      STOP 'MEMORY TEST 1'

      DO K = 1, MAIN%ISTP
!        CALL Print_SumAC2("In the time loop")
        IF (LQSTEA) THEN
          CALL QUASI_STEADY(K)
        ELSE
          CALL UN_STEADY(K,CALLFROM)
        END IF
      END DO

#ifdef TIMINGS
       CALL WAV_MY_WTIME(TIME2)
      WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----TOTAL TIME IN PROG-----', TIME2-TIME1
# endif

# if defined MPI_PARALL_GRID && !defined ROMS_WWM_PGMCL_COUPLING && !defined MODEL_COUPLING_ATM_WAV && !defined MODEL_COUPLING_OCN_WAV
      call parallel_finalize
# endif

# if defined ROMS_WWM_PGMCL_COUPLING || defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
      END SUBROUTINE
# else
      END PROGRAM
# endif
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
