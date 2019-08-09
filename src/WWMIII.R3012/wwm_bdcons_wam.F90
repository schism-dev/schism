#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_BND_INTERPOLATION_ARRAY(TheInfo)
      USE DATAPOOL
      IMPLICIT NONE
      type(FD_FORCING_GRID), intent(in) :: TheInfo
      integer IP, eIDX
      real(rkind) eX, eY
      integer eCF_IX, eCF_IY
      real(rkind) eCF_COEFF(4)
      LOGICAL EXTRAPO_OUT
      integer nbExtrapolation
      WRITE(STAT%FHNDL,*) 'Begin COMPUTE_BND_INTERPOLATION_ARRAY'
      WRITE(STAT%FHNDL,*) 'EXTRAPOLATION_ALLOWED_BOUC=', EXTRAPOLATION_ALLOWED_BOUC
      FLUSH(STAT%FHNDL)
      allocate(CF_IX_BOUC(IWBMNP), CF_IY_BOUC(IWBMNP), CF_COEFF_BOUC(4,IWBMNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('CF_*_BOUC allocation error')
      
      nbExtrapolation = 0
      DO IP=1,IWBMNP
        eIdx = IWBNDLC(IP)
        eX=XP(eIDX)
        eY=YP(eIDX)
        CALL COMPUTE_SINGLE_INTERPOLATION_INFO(TheInfo, EXTRAPOLATION_ALLOWED_BOUC, eX, eY, eCF_IX, eCF_IY, eCF_COEFF, EXTRAPO_OUT)
        CF_IX_BOUC(IP) = eCF_IX
        CF_IY_BOUC(IP) = eCF_IY
        CF_COEFF_BOUC(:,IP) = eCF_COEFF
        IF (EXTRAPO_OUT .eqv. .TRUE.) THEN
          nbExtrapolation=nbExtrapolation + 1
        END IF
      END DO
      WRITE(STAT % FHNDL, *) 'Computing extrapolation array for boundary'
      WRITE(STAT % FHNDL, *) 'nbExtrapolation=', nbExtrapolation
      END SUBROUTINE COMPUTE_BND_INTERPOLATION_ARRAY
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef GRIB_API_ECMWF
      SUBROUTINE INIT_GRIB_WAM_BOUNDARY
      USE DATAPOOL
      USE GRIB_API
      IMPLICIT NONE
      INTEGER ifile, IFILE_IN
      REAL(rkind) :: eTimeMjd
      LOGICAL STEPRANGE_IN
      LOGICAL :: USE_DATATIME = .TRUE.
      type(FD_FORCING_GRID) :: TheInfo
      character(len=20) shortName
      integer GRIB_TYPE
      integer, allocatable :: ListDir_i(:), ListFreq_i(:)
      integer nbTotalNumberEntry
      LOGICAL IsFirst
      character(len=140) eFile
      integer i, idir, ifreq, n
      integer nbdir_wam_read, nbfreq_wam_read
      integer freqScal, dirScal
      integer, allocatable :: igrib(:)
      real(rkind) eDIR, eFR, eFreq
      real(rkind) eDiff, eDiff1, eDiff2
      real(rkind) eWD1, eWD2
      logical IsAssigned
      integer IS, ID, idx
      integer ID1, ID2
      real(rkind) eTimeOut, DeltaDiff
      real(rkind) DELTH_WAM, CO1, WETAIL_WAM
      integer M
      CALL TEST_FILE_EXIST_DIE("Missing list of WAM files: ", TRIM(WAV%FNAME))
      OPEN(WAV%FHNDL,FILE=WAV%FNAME,STATUS='OLD')
      WRITE(STAT%FHNDL,*) WAV%FHNDL, WAV%FNAME, BND%FHNDL, BND%FNAME
      STEPRANGE_IN = .TRUE.
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_BOUND .or. (myrank .eq. 0)) THEN
# endif
        !
        ! Determining the number of times
        !
        NUM_WAM_SPEC_FILES = 0
        DO
          READ( WAV%FHNDL, *, IOSTAT = ISTAT )
          IF ( ISTAT /= 0 ) EXIT
          NUM_WAM_SPEC_FILES = NUM_WAM_SPEC_FILES + 1
        END DO
        REWIND(WAV%FHNDL)
        WRITE(STAT%FHNDL,*) 'NUM_WAM_SPEC_FILES=', NUM_WAM_SPEC_FILES
        !
        ! Reading the file names
        !
        ALLOCATE(WAM_SPEC_FILE_NAMES_BND(NUM_WAM_SPEC_FILES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons_wam, allocate error 9')
        DO IFILE_IN = 1, NUM_WAM_SPEC_FILES
          READ( WAV%FHNDL, *) WAM_SPEC_FILE_NAMES_BND(IFILE_IN)
        END DO
        CLOSE (WAV%FHNDL)
        !
        ! Determining total number of entries
        ! and also the number of directions and frequencies
        !
        nbTotalNumberEntry=0
        IsFirst=.TRUE.
        DO IFILE_IN = 1, NUM_WAM_SPEC_FILES
          eFile=WAM_SPEC_FILE_NAMES_BND(IFILE_IN)
          CALL TEST_FILE_EXIST_DIE("Missing grib file: ", TRIM(eFile))
          CALL GRIB_OPEN_FILE(ifile, TRIM(eFile), 'r')
          call grib_count_in_file(ifile,n)
          allocate(igrib(n))
          DO i=1,n
            call grib_new_from_file(ifile, igrib(i))
            call grib_get(igrib(i), 'numberOfDirections', nbdir_wam_read)
            call grib_get(igrib(i), 'numberOfFrequencies', nbfreq_wam_read)
            IF (IsFirst .eqv. .TRUE.) THEN
              nbdir_wam = nbdir_wam_read
              nbfreq_wam = nbfreq_wam_read
              call grib_get(igrib(i), 'directionScalingFactor', dirScal)
              call grib_get(igrib(i), 'frequencyScalingFactor', freqScal)
              allocate(ListDir_i(nbdir_wam), ListFreq_i(nbfreq_wam), ListDir_wam(nbdir_wam), ListFreq_wam(nbfreq_wam), DFIM_wam(nbFreq_wam), stat=istat)
              call grib_get(igrib(i), 'scaledDirections', ListDir_i)
              call grib_get(igrib(i), 'scaledFrequencies', ListFreq_i)
              DO idir=1,nbdir_wam
                eDir = MyREAL(ListDir_i(idir)) / MyREAL(dirScal)
                eDir = 270 - eDir
                IF (eDir .le. ZERO) THEN
                  eDir = eDir+ 360
                END IF
                IF (eDir .ge. 360) THEN
                  eDir = eDir - 360
                END IF
                ListDir_wam(idir) = eDir
                WRITE(STAT%FHNDL,*) 'idir=', idir, ' eDir=', eDir
              END DO
              DO ifreq=1,nbfreq_wam
                eFreq = MyREAL(ListFreq_i(ifreq)) / MyREAL(freqScal)
                ListFreq_wam(ifreq) = eFreq
                WRITE(STAT%FHNDL,*) 'ifreq=', ifreq, ' eFreq=', eFreq
              END DO
              FRATIO = ListFreq_wam(2) / ListFreq_wam(1)
              WRITE(STAT%FHNDL,*) 'FRATIO=', FRATIO
              DELTH_WAM = PI2 / MyREAL(nbdir_wam)
              WRITE(STAT%FHNDL,*) 'DELTH_WAM=', DELTH_WAM
              CO1 = 0.5*(FRATIO-1.)*DELTH_WAM
              WRITE(STAT%FHNDL,*) 'CO1=', CO1
              DFIM_wam(1) = CO1 * ListFreq_wam(1)
              DO M=2,nbFreq_wam-1
                 DFIM_wam(M) = CO1 * (ListFreq_wam(M) + ListFreq_wam(M-1))
              ENDDO
              DFIM_wam(nbFreq_wam) = CO1 * ListFreq_wam(nbFreq_wam-1)
              DO M=1,nbFreq_wam
                WRITE(STAT%FHNDL,*) 'M=', M, ' DFIM=', DFIM_wam(M)
              END DO
              WETAIL_WAM = 0.25
              DELT25_WAM = WETAIL_WAM*ListFreq_wam(nbFreq_wam)*DELTH_WAM
              deallocate(ListDir_i, ListFreq_i)
            ELSE
              IF ((nbdir_wam .ne. nbdir_wam_read).or.(nbfreq_wam .ne. nbfreq_wam_read)) THEN
                Print *, 'nbdir_wam =', nbdir_wam,  'nbdir_wam_read =', nbdir_wam_read
                Print *, 'nbfreq_wam=', nbfreq_wam, 'nbfreq_wam_read=', nbfreq_wam_read
                CALL WWM_ABORT('number of frequencies/directions is inconsistent')
              END IF
            END IF
            IsFirst=.FALSE.
          END DO
          deallocate(igrib)
          CALL GRIB_CLOSE_FILE(ifile)
          nbTotalNumberEntry = nbTotalNumberEntry + n / (nbdir_wam * nbfreq_wam)
        END DO
        ALLOCATE(eVAR_BOUC_WAM % ListTime(nbTotalNumberEntry), ListIFileWAM(nbTotalNumberEntry), stat=istat)
        eVAR_BOUC_WAM % nbTime = nbTotalNumberEntry
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons_wam, allocate error 9')
        idx=0
        DO IFILE_IN = 1, NUM_WAM_SPEC_FILES
          eFile=WAM_SPEC_FILE_NAMES_BND(IFILE_IN)
!          Print *, 'iFile=', iFile, ' eFile=', TRIM(eFile)
          CALL TEST_FILE_EXIST_DIE("Missing grib file: ", TRIM(eFile))
          CALL GRIB_OPEN_FILE(ifile, TRIM(eFile), 'r')
          call grib_count_in_file(ifile,n)
          allocate(igrib(n))
          DO i=1,n
            call grib_new_from_file(ifile, igrib(i))
            call grib_get(igrib(i), 'directionNumber', idir)
            call grib_get(igrib(i), 'frequencyNumber', ifreq)
            IF ((idir .eq. 1).and.(ifreq .eq. 1)) THEN
!              Print *, 'i=', i, ' idir=', idir, ' ifreq=', ifreq
              CALL RAW_READ_TIME_OF_GRIB_FILE(ifile, igrib(i), STEPRANGE_IN, eTimeOut)
              !
              idx=idx+1
              eVAR_BOUC_WAM % ListTime(idx) = eTimeOut
              ListIFileWAM(idx) = IFILE_IN
            END IF
          END DO
          deallocate(igrib)
          CALL GRIB_CLOSE_FILE(ifile)
        END DO
        !
        ! reading the grid
        !
        shortName='2dfd'
        GRIB_TYPE=1 ! 1 for ECMWF
        IFILE_IN = 1
!        Print *, 'Before READ_GRID_INFO_FROM_GRIB'
        CALL READ_GRID_INFO_FROM_GRIB(TheInfo, WAM_SPEC_FILE_NAMES_BND(IFILE_IN), shortName, GRIB_TYPE)
!        Print *, 'After READ_GRID_INFO_FROM_GRIB'
# ifdef MPI_PARALL_GRID
      END IF
# endif
      CALL COMPUTE_BND_INTERPOLATION_ARRAY(TheInfo)
      deallocate(TheInfo % LON, TheInfo % LAT)
      nx_wam = TheInfo % nx_dim
      ny_wam = TheInfo % ny_dim
!      Print *, 'After COMPUTE_BND_INTERPOLATION_ARRAY'
      !
      ! Now the spectral interpolation arrays
      !
      allocate(WAM_ID1(MDC), WAM_ID2(MDC), WAM_WD1(MDC), WAM_WD2(MDC), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('CF_*_BOUC allocation error')
      WAM_ID1=0
      WAM_ID2=0
      DO ID=1,MDC
        eDIR=SPDIR(ID) * RADDEG
        IsAssigned=.false.
        DO ID1=1,nbdir_wam
          IF (ID1 .lt. nbdir_wam) THEN
            ID2=ID1+1
          ELSE
            ID2=1
          END IF
          IF (IsAssigned .eqv. .false.) THEN
            eDiff = ListDir_wam(ID2) - ListDir_wam(ID1)
            IF (eDiff .gt. 180) THEN
              eDiff = eDiff - 360.0
            END IF
            IF (eDiff .lt. -180) THEN
              eDiff = eDiff + 360.0
            END IF
            !
            eDiff1=eDIR - ListDir_wam(ID1)
            IF (eDiff1 .gt. 180) THEN
              eDiff1 = eDiff1 - 360.0
            END IF
            IF (eDiff1 .lt. -180) THEN
              eDiff1 = eDiff1 + 360.0
            END IF
            !
            eDiff2=ListDir_wam(ID2) - eDir
            IF (eDiff2 .gt. 180) THEN
              eDiff2 = eDiff2 - 360.0
            END IF
            IF (eDiff2 .lt. -180) THEN
              eDiff2 = eDiff2 + 360.0
            END IF
            DeltaDiff = abs(eDiff) - abs(eDiff1) - abs(eDiff2)
            IF (abs(DeltaDiff) .le. 1.0) THEN
              eWD1 = eDiff2 / eDiff
              eWD2 = eDiff1 / eDiff
              IsAssigned=.TRUE.
              WAM_ID1(ID) = ID1
              WAM_ID2(ID) = ID2
              WAM_WD1(ID) = eWD1
              WAM_WD2(ID) = eWD2
            END IF
          END IF
        END DO
        IF (IsAssigned .eqv. .FALSE.) THEN
          CALL WWM_ABORT('Error in the interpolation direction')
        END IF
        WRITE(STAT%FHNDL,*) 'ID=', ID, 'eDir=', eDIR
        WRITE(STAT%FHNDL,*) 'WAM_ID12=', WAM_ID1(ID), WAM_ID2(ID)
        WRITE(STAT%FHNDL,*) 'WAM_WD12=', WAM_WD1(ID), WAM_WD2(ID)
        WRITE(STAT%FHNDL,*) 'WAM_eD12=', ListDir_wam(WAM_ID1(ID)), ListDir_wam(WAM_ID2(ID))
      END DO
      allocate(WAM_IS1(MSC), WAM_IS2(MSC), WAM_WS1(MSC), WAM_WS2(MSC), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('CF_*_BOUC allocation error')
      WAM_IS1=0
      WAM_IS2=0
      DO IS=1,MSC
        IsAssigned=.FALSE.
        eFR=FR(IS)
        DO iFreq=1,nbfreq_wam-1
          IF (IsAssigned .eqv. .FALSE.) THEN
            eDiff=ListFreq_wam(iFreq+1) - ListFreq_wam(iFreq)
            eDiff1=eFR - ListFreq_wam(iFreq)
            eDiff2=ListFreq_wam(iFreq+1) - eFR
            IF ((eDiff1 .ge. 0).and.(eDiff2 .ge.0)) THEN
              IsAssigned=.TRUE.
              WAM_IS1(IS)=iFreq
              WAM_IS2(IS)=iFreq+1
              WAM_WS1(IS)=eDiff2 / eDiff
              WAM_WS2(IS)=eDiff1 / eDiff
            END IF
          END IF
        END DO
!        WRITE(STAT%FHNDL,*) 'IS=', IS, 'eFR=', eFR
!        WRITE(STAT%FHNDL,*) 'WAM_IS12=', WAM_IS1(IS), WAM_IS2(IS)
!        WRITE(STAT%FHNDL,*) 'WAM_WS12=', WAM_WS1(IS), WAM_WS2(IS)
      END DO
!      Print *, 'Leaving INIT_GRIB_WAM_BOUNDARY'
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_GRIB_WAM_BOUNDARY_WBAC_KERNEL_NAKED(WBAC_WAM, IFILE_IN, eTimeSearch)
      USE DATAPOOL
      USE GRIB_API  
      IMPLICIT NONE
      real(rkind), intent(out) :: WBAC_WAM(nbdir_wam, nbfreq_wam, nx_wam, ny_wam)
      integer, intent(in) :: IFILE_IN
      real(rkind), intent(in) :: eTimeSearch
      !
      real(rkind) :: values(nx_wam*ny_wam)
      integer :: DirFreqStatus(nbdir_wam, nbfreq_wam)
      character(len=140) eFile
      integer i, n
      integer eDiff
      integer idx, idir, ifreq
      real(rkind) DeltaDiff
      integer, allocatable :: igrib(:)
      LOGICAL :: STEPRANGE_IN = .TRUE.
      real(rkind) eTimeOut
      character(len=140) eShortName
      integer iX, iY, ifile
!      Print *, 'Begin READ_GRIB_WAM_BOUNDARY_WBAC_KERNEL_NAKED'
      DirFreqStatus=0
      eFile=WAM_SPEC_FILE_NAMES_BND(IFILE_IN)
      CALL TEST_FILE_EXIST_DIE("Missing grib file: ", TRIM(eFile))
      CALL GRIB_OPEN_FILE(ifile, TRIM(eFile), 'r')
      call grib_count_in_file(ifile,n)
      allocate(igrib(n))
      DO i=1,n
        call grib_new_from_file(ifile, igrib(i))
        CALL RAW_READ_TIME_OF_GRIB_FILE(ifile, igrib(i), STEPRANGE_IN, eTimeOut)
        DeltaDiff = abs(eTimeOut - eTimeSearch)
        IF (DeltaDiff .le. 1.0E-8) THEN
          call grib_get(igrib(i), 'shortName', eShortName)
          IF (TRIM(eShortName) .eq. '2dfd') THEN
            call grib_get(igrib(i), 'directionNumber', idir)
            call grib_get(igrib(i), 'frequencyNumber', ifreq)
            CALL grib_get(igrib(i), 'values', values)
            DirFreqStatus(idir, ifreq) = 1
            idx=0
            DO iY=1,ny_wam
              DO iX=1,nx_wam
                idx=idx+1
                WBAC_WAM(idir, ifreq, iX,iY) = values(idx)
              END DO
            END DO  
          END IF
        END IF
      END DO
      deallocate(igrib)
      CALL GRIB_CLOSE_FILE(ifile)
      eDiff= sum(DirFreqStatus) - nbdir_wam * nbfreq_wam
      if (eDiff .ne. 0) THEN
        CALL WWM_ABORT('Error reading WAM file. Some direction/frequencies not assigned')
      END IF
      WRITE(STAT%FHNDL,*) 'sum(WBAC_WAM)=', sum(WBAC_WAM)
!      Print *, 'End READ_GRIB_WAM_BOUNDARY_WBAC_KERNEL_NAKED'
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_GRIB_WAM_BOUNDARY_WBAC_KERNEL(WBACOUT, IFILE, eTimeSearch)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(OUT)   :: WBACOUT(MSC,MDC,IWBMNP)
      integer, intent(in) :: IFILE
      real(rkind), intent(in) :: eTimeSearch
      !
      real(rkind) :: WBAC_WAM    (nbdir_wam, nbfreq_wam, nx_wam, ny_wam)
      real(rkind) :: WBAC_WAM_LOC(nbdir_wam, nbfreq_wam)
      integer ID1, ID2, IS1, IS2
      integer ID, IS, J, IP
      real(rkind) WD1, WD2, WS1, WS2
      real(rkind) ACLOC(MSC,MDC)
      integer IX, IY
      real(rkind) eAC_1, eAC_2, eAC
      real(rkind) EM, HS_WAM, eSum
      integer M, K
      LOGICAL :: DoHSchecks = .TRUE.
      real(rkind) ETOT, tmp(msc), DS, ETAIL, HS_WWM, EMwork
      
      CALL READ_GRIB_WAM_BOUNDARY_WBAC_KERNEL_NAKED(WBAC_WAM, IFILE, eTimeSearch)
      WRITE(STAT%FHNDL,*) 'RETURN: sum(WBAC_WAM)=', sum(WBAC_WAM)
      WRITE(STAT%FHNDL,*) 'IWBMNP=', IWBMNP
      DO IP=1,IWBMNP
        IX=CF_IX_BOUC(IP)
        IY=CF_IY_BOUC(IP)
        WBAC_WAM_LOC=0
        DO J=1,4
          WBAC_WAM_LOC(:,:) = WBAC_WAM_LOC(:,:) + CF_COEFF_BOUC(J,IP)*WBAC_WAM(:,:,IX+SHIFTXY(J,1),IY+SHIFTXY(J,2))
        END DO
        WRITE(STAT%FHNDL,*) 'sum(WBAC_WAM_LOC)=', sum(WBAC_WAM_LOC)
        !
        IF (DoHSchecks) THEN
          DO J=1,4
            WRITE(STAT%FHNDL,*) 'J=', J, ' eCF=', CF_COEFF_BOUC(J,IP)
          END DO
          EM=0
          DO M=1,nbfreq_wam
            eSum=0
            DO K=1,nbdir_wam
              eSum = eSum + WBAC_WAM_LOC(K,M)
            END DO
            EM = EM + DFIM_WAM(M)*eSum
            Print *, 'M=', M, ' EM=', EM
          END DO
          Print *, 'DELT25=', DELT25_WAM
          EM = EM + DELT25_WAM*eSum
          Print *, 'EM=', EM
          EMwork=MAX(ZERO, EM)
          Print *, 'EMwork=', EMwork
          HS_WAM = 4.*SQRT(EMwork)
        END IF
        ACLOC=0
        DO IS=1,MSC
          DO ID=1,MDC
            ID1=WAM_ID1(ID)
            ID2=WAM_ID2(ID)
            WD1=WAM_WD1(ID)
            WD2=WAM_WD2(ID)
            !
            IS1=WAM_IS1(ID)
            IS2=WAM_IS2(ID)
            WS1=WAM_WS1(ID)
            WS2=WAM_WS2(ID)
            !
            Print *, 'ID12=', ID1, ID2, ' IS12=', IS1, IS2
            IF (IS1 .gt. 0) THEN
              eAC_1=WD1 * WBAC_WAM_LOC(ID1, IS1) + WD2 * WBAC_WAM_LOC(ID2, IS1)
              eAC_2=WD1 * WBAC_WAM_LOC(ID1, IS2) + WD2 * WBAC_WAM_LOC(ID2, IS2)
              eAC=WS1 * eAC_1 + WS2 * eAC_2
              ACLOC(IS,ID)=eAC
            END IF
          END DO
        END DO
        IF (DoHSchecks) THEN
          ETOT=0
          DO ID=1,MDC
            tmp(:) = acloc(:,id) * spsig
            ETOT = ETOT + tmp(1) * ONEHALF * ds_incr(1)*ddir
            do is = 2, msc
              ETOT = ETOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
            end do
            ETOT = ETOT + ONEHALF * tmp(msc) * ds_incr(msc)*ddir
          END DO
          DS    = SPSIG(MSC) - SPSIG(MSC-1)
          ETAIL = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,2) * DDIR * DS
          ETOT  = ETOT + PTAIL(6) * ETAIL
          HS_WWM = 4*SQRT(MAX(0.0, ETOT))
          WRITE(STAT%FHNDL,*) 'BOUND IP=', IP, '/', IWBMNP
          WRITE(STAT%FHNDL,*) 'ETOT(WAM/WWM)=', EM, ETOT
          WRITE(STAT%FHNDL,*) 'HS(WAM/WWM)=', HS_WAM, HS_WWM, ETOT
        END IF
        WBACOUT(:,:,IP)=ACLOC
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_GRIB_WAM_BOUNDARY_WBAC(WBACOUT)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(OUT)   :: WBACOUT(MSC,MDC,IWBMNP)
      !
      integer iTime
      real(rkind) DeltaDiff, eTimeSearch
      real(rkind) eTimeDay
      integer iFile
      CHARACTER(LEN=15) :: eTimeStr
      
      eTimeSearch=MAIN % TMJD
      DO iTime=1, eVAR_BOUC_WAM % nbTime
        eTimeDay=eVAR_BOUC_WAM % ListTime(iTime)
        DeltaDiff= abs(eTimeDay - eTimeSearch)
        iFile=ListIFileWAM(iTime)
        IF (DeltaDiff .le. 1.0e-8) THEN
          CALL READ_GRIB_WAM_BOUNDARY_WBAC_KERNEL(WBACOUT, iFile, eTimeSearch)
          RETURN
        END IF
      END DO
      Print *, 'nbTime=', eVAR_BOUC_WAM % nbTime, ' eTimeSearch=', eTimeSearch
      DO iTime=1, eVAR_BOUC_WAM % nbTime
        eTimeDay=eVAR_BOUC_WAM % ListTime(iTime)
        CALL MJD2CT(eTimeDay,eTimeStr)
        Print *, 'iTime=', iTime, ' eTime=', eTimeDay, ' date=', eTimeStr
      END DO      
      CALL WWM_ABORT('Failed to find the right record in READ_GRIB_BOUNDARY_WBAC')
      END SUBROUTINE
#endif
