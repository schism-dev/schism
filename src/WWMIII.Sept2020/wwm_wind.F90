#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_WIND_INPUT
#ifdef NCDF
      USE NETCDF
#endif
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER     :: IT, IFILE
      REAL(rkind) :: WDIRT
      REAL(rkind) :: cf_w1, cf_w2

      WINDXY(:,:) = 0.0
#ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND) THEN
        MNP_WIND=MNP
        allocate(XP_WIND(MNP_WIND), YP_WIND(MNP_WIND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 1')
        XP_WIND=XP
        YP_WIND=YP
      ELSE
        MNP_WIND=np_total
        IF (myrank .eq. 0) THEN
          allocate(XP_WIND(MNP_WIND), YP_WIND(MNP_WIND), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 1')
          XP_WIND=XPtotal
          YP_WIND=YPtotal
        END IF
      END IF
#else
      MNP_WIND=MNP
      allocate(XP_WIND(MNP_WIND), YP_WIND(MNP_WIND), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 1')
      XP_WIND=XP
      YP_WIND=YP
#endif
      IF (LSTWD) THEN
        IF (LCWIN) THEN
          WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'HOMOGENOUS STEADY WIND FIELD IS USED' 
          WRITE(WINDBG%FHNDL,'("+TRACE...",A,I10)') 'WIND IS COMING FROM WWM - WINDFORMAT', IWINDFORMAT, LWDIR
          FLUSH(WINDBG%FHNDL)
          IF (LWDIR) THEN
            CALL DEG2NAUT(WDIR, WDIRT, LNAUTIN)
            WINDXY(:,1) =  WVEL * COS(WDIRT * DEGRAD)
            WINDXY(:,2) =  WVEL * SIN(WDIRT * DEGRAD)
          ELSE
            WINDXY(:,1) = CWINDX
            WINDXY(:,2) = CWINDY
          END IF
        ELSE ! LCWIN
          WRITE(WINDBG%FHNDL,'("+TRACE...",A,I10)') 'WIND IS COMING FROM WWM - WINDFORMAT', IWINDFORMAT
          WRITE(WINDBG%FHNDL,'("+TRACE...",A)')  'SPATIAL VARIABLE WIND FIELD IS USED'
          FLUSH(WINDBG%FHNDL)
          IF (IWINDFORMAT == 1) THEN
            CALL CSEVAL( WIN%FHNDL, TRIM(WIN%FNAME), .FALSE., 2, WINDXY, MULTIPLE_IN_WIND)
#ifdef NCDF
          ELSE IF (IWINDFORMAT == 2) THEN ! NETCDF created using ncl_convert2nc using DWD grib
            CALL INIT_NETCDF_DWD
            CALL FIND_WIND_NEAREST_LOWER_IDX(SEWI%BMJD, idxWind)
            IFILE=WIND_TIME_IFILE(idxWind)
            IT=WIND_TIME_IT(idxWind)
            CALL READ_NETCDF_DWD(IFILE, IT, WINDXY)
          ELSE IF (IWINDFORMAT == 3) THEN ! NETCDF created using cdo -f nc copy file.grb file.nc this is CFRS
            CALL INIT_NETCDF_CRFS
            CALL FIND_WIND_NEAREST_LOWER_IDX(SEWI%BMJD, idxWind)
            IFILE=WIND_TIME_IFILE(idxWind)
            IT=WIND_TIME_IT(idxWind)
            CALL READ_NETCDF_CRFS(IFILE, IT, WINDXY)
          ELSE IF (IWINDFORMAT == 4) THEN ! NETCDF NARR downloaded from NOMAD
            CALL INIT_NETCDF_NARR
            CALL FIND_WIND_NEAREST_LOWER_IDX(SEWI%BMJD, idxWind)
            IFILE=WIND_TIME_IFILE(idxWind)
            IT=WIND_TIME_IT(idxWind)
            CALL READ_NETCDF_NARR(IFILE, IT, WINDXY)
          ELSE IF (IWINDFORMAT == 5) THEN ! NETCDF CF_COMPLIANT STATIONARY FIELD 
            WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'COMPUTING CF INTERPOLATION COEFS AND LOADING WIND_TIME_MJD'
            FLUSH(WINDBG%FHNDL)
            CALL INIT_NETCDF_CF_WWM_WIND(eVAR_WIND)
            ALLOCATE(tmp_wind1(MNP,2),tmp_wind2(MNP,2), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 1')
            CALL GET_CF_TIME_INDEX(eVAR_WIND, REC1_wind_new,REC2_wind_new,cf_w1,cf_w2)
            CALL READ_INTERP_NETCDF_CF_WWM_WIND(REC1_wind_new,tmp_wind1)
            IF (cf_w1.NE.1) THEN
              CALL READ_INTERP_NETCDF_CF_WWM_WIND(REC2_wind_new,tmp_wind2)
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)+cf_w2*tmp_wind2(:,:)
            ELSE
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)
            END IF
          ELSE IF (IWINDFORMAT == 6) THEN ! DIRECT WWM forcing (no interp)
            CALL INIT_DIRECT_NETCDF_CF(eVAR_WIND, MULTIPLE_IN_WIND, WIN%FNAME, "Uwind")
            ALLOCATE(tmp_wind1(MNP,2),tmp_wind2(MNP,2), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 1')
            CALL GET_CF_TIME_INDEX(eVAR_WIND, REC1_wind_new,REC2_wind_new,cf_w1,cf_w2)
            CALL READ_DIRECT_NETCDF_CF(eVAR_wind, REC1_wind_new,tmp_wind1)
            IF (cf_w1.NE.1) THEN
              CALL READ_DIRECT_NETCDF_CF(eVAR_wind, REC2_wind_new,tmp_wind2)
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)+cf_w2*tmp_wind2(:,:)
            ELSE
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)
            END IF
#endif
#ifdef GRIB_API_ECMWF
          ELSE IF (IWINDFORMAT == 7) THEN ! GRIB forcing from ecmwf
            ALLOCATE(tmp_wind1(MNP,2),tmp_wind2(MNP,2), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 1')
            CALL GET_CF_TIME_INDEX(eVAR_WIND, REC1_wind_new,REC2_wind_new,cf_w1,cf_w2)
            CALL READ_GRIB_WIND(REC1_wind_new,tmp_wind1)
            IF (cf_w1.NE.1) THEN
              CALL READ_GRIB_WIND(REC2_wind_new,tmp_wind2)
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)+cf_w2*tmp_wind2(:,:)
            ELSE
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)
            END IF
#endif
          ELSE
            CALL wwm_abort('Wrong choice of IWINDFORMAT (maybe need NETCDF or GRIB)')
          END IF
        ENDIF
      ELSE IF (LSEWD) THEN
        IF (LCWIN) THEN
          CALL wwm_abort('LSEWD + LCWIN NOT READY')
!         CALL READ_WIND_TIME_SERIES(IT) ! set time according to wwminput.nml and get initial time step
!         CALL SET_INITIAL_WIND(IT) ! 
        ELSE
          WRITE(WINDBG%FHNDL,'("+TRACE...",A,I10)') 'WIND IS COMING FROM WWM - WINDFORMAT', IWINDFORMAT
          WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'NONSTATIONARY WIND FIELD IS USED        '
          FLUSH(WINDBG%FHNDL)
          SEWI%TOTL = (SEWI%EMJD - SEWI%BMJD) * DAY2SEC
          SEWI%ISTP = NINT( SEWI%TOTL / SEWI%DELT ) + 1
          SEWI%TMJD = SEWI%BMJD
          WRITE(WINDBG%FHNDL,*) SEWI%BEGT, SEWI%ENDT, SEWI%ISTP, SEWI%TOTL/3600.0, SEWI%DELT
          WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'SPATIAL VARIABLE WIND FIELD IS USED'
          WRITE(WINDBG%FHNDL,*) 'IWINDFORMAT=', IWINDFORMAT
          FLUSH(WINDBG%FHNDL)
          IF (IWINDFORMAT == 1) THEN
            OPEN(WIN%FHNDL, FILE = TRIM(WIN%FNAME), STATUS = 'OLD')
            CALL CSEVAL( WIN%FHNDL, TRIM(WIN%FNAME), .TRUE., 2, WINDXY, MULTIPLE_IN_WIND)
#ifdef NCDF
          ELSE IF (IWINDFORMAT == 2) THEN ! NETCDF created using ncl_convert2nc using DWD grib
            CALL INIT_NETCDF_DWD
            CALL FIND_WIND_NEAREST_LOWER_IDX(MAIN%TMJD, idxWind)
            IFILE=WIND_TIME_IFILE(idxWind)
            IT=WIND_TIME_IT(idxWind)
            CALL READ_NETCDF_DWD(IFILE, IT, WINDXY)
          ELSE IF (IWINDFORMAT == 3) THEN ! NETCDF created using cdo -f nc copy file.grb file.nc
            CALL INIT_NETCDF_CRFS
            CALL FIND_WIND_NEAREST_LOWER_IDX(MAIN%TMJD, idxWind)
            IFILE=WIND_TIME_IFILE(idxWind)
            IT=WIND_TIME_IT(idxWind)
            CALL READ_NETCDF_CRFS(IFILE, IT, WINDXY)
          ELSE IF (IWINDFORMAT == 4) THEN ! NETCDF created using cdo -f nc copy file.grb file.nc
            CALL INIT_NETCDF_NARR
            CALL FIND_WIND_NEAREST_LOWER_IDX(MAIN%TMJD, idxWind)
            IFILE=WIND_TIME_IFILE(idxWind)
            IT=WIND_TIME_IT(idxWind)
            CALL READ_NETCDF_NARR(IFILE, IT, WINDXY)
          ELSE IF (IWINDFORMAT == 5) THEN
            WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'SPATIAL/TEMPORAL VARIABLE WIND FIELD IS USED CF NETCDF'
            WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'COMPUTING CF INTERPOLATION COEFS AND LOADING WIND_TIME_MJD'
            FLUSH(WINDBG%FHNDL)
            CALL INIT_NETCDF_CF_WWM_WIND(eVAR_WIND)
            ALLOCATE(tmp_wind1(MNP,2), tmp_wind2(MNP,2), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 2')
            CALL GET_CF_TIME_INDEX(eVAR_WIND, REC1_wind_new,REC2_wind_new,cf_w1,cf_w2)
            CALL READ_INTERP_NETCDF_CF_WWM_WIND(REC1_wind_new,tmp_wind1)
            IF (cf_w1.NE.1) THEN
              CALL READ_INTERP_NETCDF_CF_WWM_WIND(REC2_wind_new,tmp_wind2)
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)+cf_w2*tmp_wind2(:,:)
            ELSE
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)
            END IF
          ELSE IF (IWINDFORMAT == 6) THEN
            WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'SPATIAL/TEMPORAL VARIABLE WIND FIELD IS USED CF NETCDF'
            WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'LOADING WIND_TIME_MJD DEFINED AT NODES'
            FLUSH(WINDBG%FHNDL)
            CALL INIT_DIRECT_NETCDF_CF(eVAR_WIND, MULTIPLE_IN_WIND, WIN%FNAME, "Uwind")
            ALLOCATE(tmp_wind1(MNP,2), tmp_wind2(MNP,2), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 2')
            CALL GET_CF_TIME_INDEX(eVAR_WIND, REC1_wind_new,REC2_wind_new,cf_w1,cf_w2)
            CALL READ_DIRECT_NETCDF_CF(eVAR_wind, REC1_wind_new,tmp_wind1)
            IF (cf_w1.NE.1) THEN
              CALL READ_DIRECT_NETCDF_CF(eVAR_wind, REC2_wind_new,tmp_wind2)
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)+cf_w2*tmp_wind2(:,:)
            ELSE
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)
            END IF
#endif
#ifdef GRIB_API_ECMWF
          ELSE IF (IWINDFORMAT == 7) THEN
            CALL INIT_GRIB_WIND
            ALLOCATE(tmp_wind1(MNP,2), tmp_wind2(MNP,2), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 2')
            CALL GET_CF_TIME_INDEX(eVAR_WIND, REC1_wind_new,REC2_wind_new,cf_w1,cf_w2)
            CALL READ_GRIB_WIND(REC1_wind_new,tmp_wind1)
            IF (cf_w1.NE.1) THEN
              CALL READ_GRIB_WIND(REC2_wind_new,tmp_wind2)
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)+cf_w2*tmp_wind2(:,:)
            ELSE
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)
            END IF
#endif
          ENDIF
        ENDIF
      ENDIF
      write(WINDBG%FHNDL,'("+TRACE... Done with CF init, Uwind ",F7.2,2x,F7.2)')minval(WINDXY(:,1)),maxval(WINDXY(:,1))
      write(WINDBG%FHNDL,'("+TRACE... Done with CF init, Vwind ",F7.2,2x,F7.2)')minval(WINDXY(:,2)),maxval(WINDXY(:,2))
      FLUSH(WINDBG%FHNDL)
 
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE UPDATE_WIND(K)
#ifdef NCDF
      USE NETCDF 
#endif
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind)             :: TMP(MNP,2)
#if defined NCDF || defined GRIB_API_ECMWF
      REAL(rkind)             :: cf_w1, cf_w2
      INTEGER                 :: IT, IFILE
#endif
      INTEGER, intent(in)     :: K
!AR: All crap ... defining K without using means that nobody has ever checked the results or anything else, so why coding at all?
!AR: Mathieu can you please fix this !!!

      WRITE(WINDBG%FHNDL,*) 'MAIN%TMJD=', MAIN%TMJD
      WRITE(WINDBG%FHNDL,*) 'SEWI(TMJD,EMJD)=', SEWI%TMJD, SEWI%EMJD
      IF ( LSEWD .AND. (MAIN%TMJD .ge. SEWI%TMJD-1.E-8) .AND. (MAIN%TMJD .le. SEWI%EMJD+1.e-8) ) THEN
        IF (IWINDFORMAT == 1) THEN
!NDM: Need to add the facility for LINTERWD
          CALL CSEVAL( WIN%FHNDL, WIN%FNAME, .TRUE., 2, TMP, MULTIPLE_IN_WIND)
          DVWIND = (TMP-WINDXY)/SEWI%DELT*MAIN%DELT
#ifdef NCDF
        ELSE IF (IWINDFORMAT == 2) THEN ! DWD_NETCDF
          CALL MOVE_BY_ONE_INDEX(IFILE, IT)
          CALL READ_NETCDF_DWD(IFILE, IT, TMP)
          DVWIND = (TMP-WINDXY)/SEWI%DELT*MAIN%DELT
        ELSE IF (IWINDFORMAT == 3) THEN ! NOAA CFRS ... the 1st step is analysis and then we have 5 + 1 forecasts, which give one the option to use either only the 6 forecast's after the analysis or use the analysis with 5 forecast's
          CALL MOVE_BY_ONE_INDEX(IFILE, IT)
          CALL READ_NETCDF_CRFS(IFILE, IT, TMP)
          DVWIND = (TMP-WINDXY)/SEWI%DELT*MAIN%DELT
        ELSE IF (IWINDFORMAT == 4) THEN ! NOAA NARR ...
          CALL MOVE_BY_ONE_INDEX(IFILE, IT)
          CALL READ_NETCDF_NARR(IFILE, IT, TMP)
          DVWIND = (TMP-WINDXY)/SEWI%DELT*MAIN%DELT
        ELSE IF (IWINDFORMAT == 5) THEN
          IF (K.EQ.1) THEN
            REC1_wind_old = 0
            REC2_wind_old = 0
          END IF
          CALL GET_CF_TIME_INDEX(eVAR_WIND, REC1_wind_new,REC2_wind_new,cf_w1,cf_w2)
          IF (REC1_wind_new.NE.REC1_wind_old) THEN
            CALL READ_INTERP_NETCDF_CF_WWM_WIND(REC1_wind_new,tmp_wind1)
          END IF
          IF (REC2_wind_new.NE.REC2_wind_old) THEN
            CALL READ_INTERP_NETCDF_CF_WWM_WIND(REC2_wind_new,tmp_wind2)
          END IF
          IF (cf_w1.NE.1) THEN
            WINDXY(:,:) = cf_w1*tmp_wind1(:,:)+cf_w2*tmp_wind2(:,:)
          ELSE
            WINDXY(:,:) = cf_w1*tmp_wind1(:,:)
          END IF
          REC1_wind_old = REC1_wind_new
          REC2_wind_old = REC2_wind_new
        ELSE IF (IWINDFORMAT == 6) THEN
          IF (K.EQ.1) THEN
            REC1_wind_old = 0
            REC2_wind_old = 0
          END IF
          CALL GET_CF_TIME_INDEX(eVAR_WIND, REC1_wind_new,REC2_wind_new,cf_w1,cf_w2)
          IF (REC1_wind_new.NE.REC1_wind_old) THEN
            CALL READ_DIRECT_NETCDF_CF(eVAR_wind, REC1_wind_new,tmp_wind1)
          END IF
          IF (REC2_wind_new.NE.REC2_wind_old) THEN
            CALL READ_DIRECT_NETCDF_CF(eVAR_wind, REC2_wind_new,tmp_wind2)
          END IF
          IF (cf_w1.NE.1) THEN
            WINDXY(:,:) = cf_w1*tmp_wind1(:,:)+cf_w2*tmp_wind2(:,:)
          ELSE
            WINDXY(:,:) = cf_w1*tmp_wind1(:,:)
          END IF
          REC1_wind_old = REC1_wind_new
          REC2_wind_old = REC2_wind_new
#endif
#ifdef GRIB_API_ECMWF
        ELSE IF (IWINDFORMAT == 7) THEN
          IF (K.EQ.1) THEN
            REC1_wind_old = 0
            REC2_wind_old = 0
          END IF
          CALL GET_CF_TIME_INDEX(eVAR_WIND, REC1_wind_new,REC2_wind_new,cf_w1,cf_w2)
          IF (REC1_wind_new.NE.REC1_wind_old) THEN
            CALL READ_GRIB_WIND(REC1_wind_new,tmp_wind1)
          END IF
          IF (REC2_wind_new.NE.REC2_wind_old) THEN
            CALL READ_GRIB_WIND(REC2_wind_new,tmp_wind2)
          END IF
          IF (cf_w1.NE.1) THEN
            WINDXY(:,:) = cf_w1*tmp_wind1(:,:)+cf_w2*tmp_wind2(:,:)
          ELSE
            WINDXY(:,:) = cf_w1*tmp_wind1(:,:)
          END IF
          REC1_wind_old = REC1_wind_new
          REC2_wind_old = REC2_wind_new
#endif
        END IF
        SEWI%TMJD = SEWI%TMJD + SEWI%DELT*SEC2DAY
      END IF
      write(WINDBG%FHNDL,'("max WINDXY:",2F7.2)')maxval(WINDXY(:,1)),maxval(WINDXY(:,2))
      write(WINDBG%FHNDL,'("min WINDXY:",2F7.2)')minval(WINDXY(:,1)),minval(WINDXY(:,2))

      IF (LWINDSWAN) THEN
        WRITE(3333,*) SEWI%TMJD
        WRITE(3333,*) WINDXY(:,1)
        WRITE(3333,*) WINDXY(:,2)
      END IF
      IF (LSEWD.AND.(IWINDFORMAT.NE.5).AND.(IWINDFORMAT.NE.6) ) THEN
        WINDXY = WINDXY + DVWIND
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MOVE_BY_ONE_INDEX(IFILE, IT)
      USE DATAPOOL
      implicit none
      integer, intent(out) :: IFILE, IT
      idxWind =idxWind+1
      IF (idxWind .gt. NDT_WIND_ALL_FILES) THEN
        CALL WWM_ABORT('Need wind after the time')
      END IF
      IFILE=WIND_TIME_IFILE(idxWind)
      IT=WIND_TIME_IT(idxWind)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE FIND_WIND_NEAREST_LOWER_IDX(eTime, idx)
      USE DATAPOOL
      implicit none
      real(rkind), intent(in) :: eTime
      integer, intent(out) :: idx
      CHARACTER(LEN=15) :: eTimeStr
      integer eIdxF, eIdx
      eIdxF=-1
      DO eIdx=1,NDT_WIND_ALL_FILES
        IF (WIND_TIME_ALL_FILES(eIdx) .le. eTime + THR8) THEN
          eIdxF=eIdx
        ENDIF
      END DO
      IF (eIdxF .eq. -1) THEN
        WRITE(WINDBG%FHNDL,*) 'NDT_WIND_ALL_FILES=', NDT_WIND_ALL_FILES
        DO eIdx=1,NDT_WIND_ALL_FILES
          CALL MJD2CT(WIND_TIME_ALL_FILES(eIdx),eTimeStr)
          WRITE(WINDBG%FHNDL,*) ' eIdx=', eIdx
          WRITE(WINDBG%FHNDL,*) ' eTime=', WIND_TIME_ALL_FILES(eIdx)
          WRITE(WINDBG%FHNDL,*) ' eTimeStr=', eTimeStr
        END DO
        CALL FLUSH(WINDBG%FHNDL)
        CALL WWM_ABORT('We failed to find the wind index')
      END IF
      idx=eIdxF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE KERNEL_INTERP_UV_WINDFD(outwind)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER I, J
      REAL(rkind), INTENT(out)           :: outwind(MNP_WIND,2)
      REAL(rkind) :: Uw, Vw
      INTEGER IX, IY
      REAL(rkind) :: cf_scale_factor, cf_add_offset
      cf_scale_factor = eVAR_WIND % cf_scale_factor
      cf_add_offset = eVAR_WIND % cf_add_offset
      DO I = 1, MNP_WIND
        Uw=ZERO
        Vw=ZERO
        IX=CF_IX(I)
        IY=CF_IY(I)
        DO J=1,4
          Uw=Uw + CF_COEFF(J,I)*UWIND_FD(IX+SHIFTXY(J,1),IY+SHIFTXY(J,2))
          Vw=Vw + CF_COEFF(J,I)*VWIND_FD(IX+SHIFTXY(J,1),IY+SHIFTXY(J,2))
        END DO
        outwind(I,1)=Uw*cf_scale_factor + cf_add_offset
        outwind(I,2)=Vw*cf_scale_factor + cf_add_offset
      END DO
      WRITE(WINDBG%FHNDL,*) 'KERNEL_INTERP_UV_WINDFD'
      WRITE(WINDBG%FHNDL,*) 'UWIND_FD, min/max=', minval(UWIND_FD), maxval(UWIND_FD)
      WRITE(WINDBG%FHNDL,*) 'VWIND_FD, min/max=', minval(VWIND_FD), maxval(VWIND_FD)
      WRITE(WINDBG%FHNDL,*) 'UWIND_FE, min/max=', minval(outwind(:,1)), maxval(outwind(:,1))
      WRITE(WINDBG%FHNDL,*) 'VWIND_FE, min/max=', minval(outwind(:,2)), maxval(outwind(:,2))
!      WRITE(WINDBG%FHNDL,*) 'max(CF_COEFF)=', maxval(abs(CF_COEFF))
!      WRITE(WINDBG%FHNDL,*) 'cf_scale_factor=', cf_scale_factor
!      WRITE(WINDBG%FHNDL,*) 'cf_add_offset=', cf_add_offset

      FLUSH(WINDBG%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE LOAD_INTERP_ARRAY(FileSave, success)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
      logical, intent(out) :: success
      character(len=256), intent(in) :: FileSave
      character (len = *), parameter :: CallFct = "LOAD_INTERP_ARRAY"
      integer, allocatable :: CF_IX_GLOBAL(:), CF_IY_GLOBAL(:)
      real(rkind), allocatable :: CF_COEFF_GLOBAL(:,:)
      integer, allocatable :: ListFirstMNP(:)
      integer, allocatable :: CF_IX_loc(:), CF_IY_loc(:)
      real(rkind), allocatable :: CF_COEFF_loc(:,:)
      integer iret, ncid, varid
      integer IP, IPglob, iPROC, NPloc, IPloc
      INQUIRE(FILE=TRIM(FileSave), EXIST=LPRECOMP_EXIST)
      IF (LPRECOMP_EXIST .eqv. .FALSE.) THEN
        success=.FALSE.
        RETURN
      END IF
      success=.TRUE.
#ifdef MPI_PARALL_GRID
      IF (myrank .eq. 0) THEN
#endif
        allocate(CF_IX_GLOBAL(np_total), CF_IY_GLOBAL(np_total), CF_COEFF_GLOBAL(4,np_total), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 52')
        !
        iret=nf90_open(TRIM(FileSave), NF90_NOWRITE, ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
        !
        iret=nf90_inq_varid(ncid, "CF_IX", varid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)
        iret=NF90_GET_VAR(ncid, varid, CF_IX_GLOBAL)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)
        !
        iret=nf90_inq_varid(ncid, "CF_IY", varid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)
        iret=NF90_GET_VAR(ncid, varid, CF_IY_GLOBAL)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)
        !
        iret=nf90_inq_varid(ncid, "CF_COEFF", varid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)
        iret=NF90_GET_VAR(ncid, varid, CF_COEFF_GLOBAL)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)
        !      
        iret=nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 27, iret)
        !
#ifdef MPI_PARALL_GRID
      END IF
#endif
      !
#ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .eqv. .FALSE.) THEN
        CF_IX=CF_IX_GLOBAL
        CF_IY=CF_IY_GLOBAL
        CF_COEFF=CF_COEFF_GLOBAL
        deallocate(CF_IX_GLOBAL, CF_IY_GLOBAL, CF_COEFF_GLOBAL)
      ELSE
        IF (myrank .eq. 0) THEN
          allocate(ListFirstMNP(nproc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 52')
          ListFirstMNP=0
          DO iProc=2,nproc
            ListFirstMNP(iProc)=ListFirstMNP(iProc-1) + ListMNP(iProc-1)
          END DO
          DO IP=1,MNP
            IPglob=iplg(IP)
            CF_IX(IP)=CF_IX_GLOBAL(IPglob)
            CF_IY(IP)=CF_IY_GLOBAL(IPglob)
            CF_COEFF(:,IP)=CF_COEFF_GLOBAL(:,IPglob)
          END DO
          !
          DO iPROC=2,nproc
            NPloc=ListMNP(iPROC)
            allocate(CF_IX_loc(NPloc), CF_IY_loc(NPloc), CF_COEFF_loc(4,NPloc), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 52')
            DO IPloc=1,NPloc
              IPglob=ListIPLG(IPloc + ListFirstMNP(iPROC))
              CF_IX_loc(IPloc)=CF_IX_GLOBAL(IPglob)
              CF_IY_loc(IPloc)=CF_IY_GLOBAL(IPglob)
              CF_COEFF_loc(:, IPloc)=CF_COEFF_GLOBAL(:, IPglob)
            END DO
            CALL MPI_SEND(CF_IX_loc, NPloc, itype, iPROC-1, 711, comm, ierr)
            CALL MPI_SEND(CF_IY_loc, NPloc, itype, iPROC-1, 712, comm, ierr)
            CALL MPI_SEND(CF_COEFF_loc, 4*NPloc, rtype, iPROC-1, 713, comm, ierr)
            deallocate(CF_IX_loc, CF_IY_loc, CF_COEFF_loc)
          END DO
          deallocate(ListFirstMNP)
          deallocate(CF_IX_GLOBAL, CF_IY_GLOBAL, CF_COEFF_GLOBAL)
        ELSE
          CALL MPI_RECV(CF_IX, MNP, itype, 0, 711, comm, istatus, ierr)
          CALL MPI_RECV(CF_IY, MNP, itype, 0, 712, comm, istatus, ierr)
          CALL MPI_RECV(CF_COEFF, 4*MNP, rtype, 0, 713, comm, istatus, ierr)
        END IF
      END IF
#else
      CF_IX=CF_IX_GLOGAL
      CF_IY=CF_IY_GLOBAL
      CF_COEFF=CF_COEFF_GLOBAL
      deallocate(CF_IX_GLOBAL, CF_IY_GLOBAL, CF_COEFF_GLOBAL)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SAVE_INTERP_ARRAY(FileSave)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
      character(len=256), intent(in) :: FileSave
      character (len = *), parameter :: CallFct = "SAVE_INTERP_ARRAY"
      integer, allocatable :: CF_IX_GLOBAL(:), CF_IY_GLOBAL(:)
      real(rkind), allocatable :: CF_COEFF_GLOBAL(:,:)
      integer, allocatable :: CF_IX_loc(:), CF_IY_loc(:)
      real(rkind), allocatable :: CF_COEFF_loc(:,:)
      integer, allocatable :: ListFirstMNP(:)
      integer ncid, iret, var_id
      integer mnp_dims, four_dims
      integer IP, IPglob, iPROC, NP_RESloc, IPloc
      WRITE(STAT%FHNDL,*) 'minval(CF_IX)=', minval(CF_IX)
      WRITE(STAT%FHNDL,*) 'minval(CF_IY)=', minval(CF_IY)
      WRITE(STAT%FHNDL,*) 'minval(CF_COEFF)=', minval(CF_COEFF)
#ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .eqv. .TRUE.) THEN
        IF (myrank .eq. 0) THEN
          allocate(ListFirstMNP(nproc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 52')
          ListFirstMNP=0
          DO iProc=2,nproc
            ListFirstMNP(iProc)=ListFirstMNP(iProc-1) + ListMNP(iProc-1)
          END DO
          !
          allocate(CF_IX_GLOBAL(np_total), CF_IY_GLOBAL(np_total), CF_COEFF_GLOBAL(4,np_total), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 52')
          !
          DO IP=1,NP_RES
            IPglob=iplg(IP)
            CF_IX_GLOBAL(IPglob)=CF_IX(IP)
            CF_IY_GLOBAL(IPglob)=CF_IY(IP)
            CF_COEFF_GLOBAL(:, IPglob)=CF_COEFF(:,IP)
          END DO
          DO iPROC=2,nproc
            NP_RESloc=ListNP_RES(iPROC)
            WRITE(STAT%FHNDL,*) ' iPROC=', iPROC, ' NP_RES_loc=', NP_RESloc 
            allocate(CF_IX_loc(NP_RESloc), CF_IY_loc(NP_RESloc), CF_COEFF_loc(4,NP_RESloc), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 52')
            WRITE(STAT%FHNDL,*) 'Step 1'
            CALL MPI_RECV(CF_IX_loc, NP_RESloc, itype, iProc-1, 611, comm, istatus, ierr)
            WRITE(STAT%FHNDL,*) 'Step 2'
            CALL MPI_RECV(CF_IY_loc, NP_RESloc, itype, iProc-1, 612, comm, istatus, ierr)
            WRITE(STAT%FHNDL,*) 'Step 3'
            CALL MPI_RECV(CF_COEFF_loc, 4*NP_RESloc, rtype, iProc-1, 613, comm, istatus, ierr)
            WRITE(STAT%FHNDL,*) 'Step 4'
            DO IPloc=1,NP_RESloc
              IPglob=ListIPLG(IPloc + ListFirstMNP(iProc))
              CF_IX_GLOBAL(IPglob)=CF_IX_loc(IPloc)
              CF_IY_GLOBAL(IPglob)=CF_IY_loc(IPloc)
              CF_COEFF_GLOBAL(:, IPglob)=CF_COEFF_loc(:, IPloc)
            END DO
            deallocate(CF_IX_loc, CF_IY_loc, CF_COEFF_loc)
          END DO
          deallocate(ListFirstMNP)
        ELSE
          WRITE(STAT%FHNDL,*) ' NP_RES=', NP_RES
          CALL MPI_SEND(CF_IX, NP_RES, itype, 0, 611, comm, ierr)
          WRITE(STAT%FHNDL,*) 'Step 1'
          CALL MPI_SEND(CF_IY, NP_RES, itype, 0, 612, comm, ierr)
          WRITE(STAT%FHNDL,*) 'Step 2'
          CALL MPI_SEND(CF_COEFF, 4*NP_RES, rtype, 0, 613, comm, ierr)
          WRITE(STAT%FHNDL,*) 'Step 3'
        END IF
      ELSE
        allocate(CF_IX_GLOBAL(np_total), CF_IY_GLOBAL(np_total), CF_COEFF_GLOBAL(4,np_total), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 52')
        CF_IX_GLOBAL=CF_IX
        CF_IY_GLOBAL=CF_IY
        CF_COEFF_GLOBAL=CF_COEFF
      END IF
#else
      allocate(CF_IX_GLOBAL(np_total), CF_IY_GLOBAL(np_total), CF_COEFF_GLOBAL(4,np_total), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 52')
      CF_IX_GLOBAL=CF_IX
      CF_IY_GLOBAL=CF_IY
      CF_COEFF_GLOBAL=CF_COEFF
#endif
      !
      ! Now writing up
      !
#ifdef MPI_PARALL_GRID
      IF (myrank .eq. 0) THEN
#endif
        iret=nf90_create(TRIM(FileSave), nf90_CLOBBER, ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
        !
        iret = nf90_def_dim(ncid, 'mnp', np_total, mnp_dims)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)
        !
        iret = nf90_def_dim(ncid, 'four', 4, four_dims)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)
        !
        iret=nf90_def_var(ncid,'CF_IX',NF90_INT,(/ mnp_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
        !
        iret=nf90_def_var(ncid,'CF_IY',NF90_INT,(/ mnp_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
        !
        iret=nf90_def_var(ncid,'CF_COEFF',NF90_RUNTYPE,(/ four_dims, mnp_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
        !
        iret=nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, iret)
        !
        ! Now writing the data
        !
        WRITE(STAT%FHNDL,*) 'minval(CF_IX_GLOBAL)=', minval(CF_IX_GLOBAL)
        WRITE(STAT%FHNDL,*) 'minval(CF_IY_GLOBAL)=', minval(CF_IY_GLOBAL)
        WRITE(STAT%FHNDL,*) 'minval(CF_COEFF_GLOBAL)=', minval(CF_COEFF_GLOBAL)
        !
        iret=nf90_open(TRIM(FileSave), NF90_WRITE, ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, iret)
        !
        iret=nf90_inq_varid(ncid,'CF_IX',var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, iret)
        !
        iret=nf90_put_var(ncid,var_id,CF_IX_GLOBAL,start = (/ 1 /), count=(/np_total/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, iret)
        !
        iret=nf90_inq_varid(ncid,'CF_IY',var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, iret)
        !
        iret=nf90_put_var(ncid,var_id,CF_IY_GLOBAL,start = (/ 1 /), count=(/np_total/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, iret)
        !
        iret=nf90_inq_varid(ncid,'CF_COEFF',var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, iret)
        !
        iret=nf90_put_var(ncid,var_id,CF_COEFF_GLOBAL,start = (/ 1, 1 /), count=(/4, np_total/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, iret)
        !
        iret=nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 15, iret)
        !
        deallocate(CF_IX_GLOBAL, CF_IY_GLOBAL, CF_COEFF_GLOBAL)
#ifdef MPI_PARALL_GRID
      END IF
#endif      
         
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_SINGLE_INTERPOLATION_INFO(TheInfo, EXTRAPO_IN, eX, eY, eCF_IX, eCF_IY, eCF_COEFF, EXTRAPO_OUT)
      USE DATAPOOL
      IMPLICIT NONE
      type(FD_FORCING_GRID), intent(in) :: TheInfo
      logical, intent(in) :: EXTRAPO_IN
      real(rkind), intent(in) :: eX, eY
      integer, intent(out) :: eCF_IX, eCF_IY
      real(rkind), intent(out) :: eCF_COEFF(4)
      logical, intent(out) :: EXTRAPO_OUT
      !
      integer IX, IY
      integer IXs, IYs
      integer IXmin, IYmin, IXmax, IYmax
      integer nx, ny
      integer aShift
      REAL(rkind) :: WI(3), X(3), Y(3), a, b
      REAL(rkind) :: MinDist, eDist
      nx = TheInfo % nx_dim
      ny = TheInfo % ny_dim
      MinDist=LARGE
!      WRITE(WINDBG%FHNDL,*) 'Start finding closest coordinate'
!      FLUSH(WINDBG%FHNDL)
      DO IX=1,nx-1
        DO IY=1,ny-1
          eDist=(eX-TheInfo % LON(IX,IY))**2 + (eY-TheInfo % LAT(IX,IY))**2
          IF (eDist .lt. MinDist) THEN
            MinDist=eDist
            IXs=IX
            IYs=IY
          END IF
        END DO
      END DO
      aShift=1
!      WRITE(WINDBG%FHNDL,*) 'lon(eX)=', eX
!      WRITE(WINDBG%FHNDL,*) 'lat(eY)=', eY
!      WRITE(WINDBG%FHNDL,*) 'LON(IXs+-1,IYs)=', TheInfo % LON(IXs-1,IYs), TheInfo % LON(IXs,IYs),TheInfo % LON(IXs+1,IYs)
!      WRITE(WINDBG%FHNDL,*) 'LAT(IXs,IYs+-1)=', TheInfo % LAT(IXs,IYs-1), TheInfo % LAT(IXs,IYs),TheInfo % LAT(IXs,IYs+1)
!      WRITE(WINDBG%FHNDL,*) 'MinDist(km)=', SQRT(MinDist)*110
!      FLUSH(WINDBG%FHNDL)

      EXTRAPO_OUT=.TRUE.
      DO
        IXmin=max(1, IXs - aShift)
        IYmin=max(1, IYs - aShift)
        IXmax=min(nx-1, IXs+aShift)
        IYmax=min(ny-1, IYs+aShift)
        DO IX=IXmin,IXmax
          DO IY=IYmin,IYmax
            ! 
            ! First triangle
            ! 
            X(1)=TheInfo % LON(IX, IY)
            X(2)=TheInfo % LON(IX+1, IY)
            X(3)=TheInfo % LON(IX, IY+1)
            Y(1)=TheInfo % LAT(IX, IY)
            Y(2)=TheInfo % LAT(IX+1, IY)
            Y(3)=TheInfo % LAT(IX, IY+1)
            CALL INTELEMENT_COEF(X,Y,eX,eY,WI)
            IF (minval(WI) .ge. -THR) THEN
              EXTRAPO_OUT=.FALSE.
              eCF_IX=IX
              eCF_IY=IY
              a=WI(2)
              b=WI(3)
              eCF_COEFF(1)=(1-a)*(1-b)
              eCF_COEFF(2)=a*(1-b)
              eCF_COEFF(3)=(1-a)*b
              eCF_COEFF(4)=a*b
              RETURN
            END IF
            !
            ! Second triangle
            !
            X(1)=TheInfo % LON(IX+1, IY+1)
            X(2)=TheInfo % LON(IX+1, IY)
            X(3)=TheInfo % LON(IX, IY+1)
            Y(1)=TheInfo % LAT(IX+1, IY+1)
            Y(2)=TheInfo % LAT(IX+1, IY)
            Y(3)=TheInfo % LAT(IX, IY+1)
            CALL INTELEMENT_COEF(X,Y,eX,eY,WI)
            IF (minval(WI) .ge. -THR) THEN
              EXTRAPO_OUT=.FALSE.
              eCF_IX=IX
              eCF_IY=IY
              a=1 - WI(3)
              b=1 - WI(2)
              eCF_COEFF(1)=(1-a)*(1-b)
              eCF_COEFF(2)=a*(1-b)
              eCF_COEFF(3)=(1-a)*b
              eCF_COEFF(4)=a*b
            END IF
          END DO
        END DO
        IF ((IXmin .eq. 1).and.(IYmin .eq. 1).and.(IXmax .eq. nx-1).and.(IYmax .eq. ny-1)) THEN
          EXIT
        END IF
        aShift=aShift + 1
      END DO
      IF ( (EXTRAPO_IN .eqv. .FALSE.) .and. (EXTRAPO_OUT .eqv. .TRUE.) ) THEN
        WRITE(STAT % FHNDL,*) 'aShift=', aShift
        WRITE(STAT % FHNDL,*) 'eX=', eX, 'eY=', eY
        FLUSH(STAT % FHNDL)
        CALL WWM_ABORT('We find a model point outside of the available forcing grid')
      ELSE
        eCF_IX = IXs
        eCF_IY = IYs
        eCF_COEFF(1)=1
        eCF_COEFF(2)=0
        eCF_COEFF(3)=0
        eCF_COEFF(4)=0
!        EXTRAPO_OUT=.TRUE.
!        WRITE(STAT % FHNDL,*) 'Point ', eX, '/', eY, ' outside grid'
!        WRITE(STAT % FHNDL,*) 'MinDist=', MinDist
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_CF_COEFFICIENTS(TheInfo)
      USE DATAPOOL
      IMPLICIT NONE
      type(FD_FORCING_GRID), intent(in) :: TheInfo
      integer I, IX, IY, IXs, IYs, IXmin, IYmin, IXmax, IYmax
      integer aShift, WeFind
      real(rkind) eDist, MinDist
      real(rkind), allocatable :: dist(:,:)
      real(rkind) closest_r(2)
      integer     closest(2)
      real(rkind) d_lon, d_lat
      integer i11, j11, i12, j12, i21, j21
      integer eCF_IX, eCF_IY
      real(rkind) eCF_COEFF(4)
      integer :: nbExtrapolation = 0
      real(rkind) :: MaxMinDist = 0
      character(len=256) :: FileSave = "wwm_filesave_interp_array.nc"
      logical success
      logical EXTRAPO_OUT
      real(rkind) eX, eY
      WRITE(WINDBG%FHNDL,*) 'Starting node loop for calcs of coefs'
      allocate(CF_IX(MNP_WIND), CF_IY(MNP_WIND), SHIFTXY(4,2), CF_COEFF(4,MNP_WIND), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 52')
      SHIFTXY(1,1)=0
      SHIFTXY(1,2)=0
      SHIFTXY(2,1)=1
      SHIFTXY(2,2)=0
      SHIFTXY(3,1)=0
      SHIFTXY(3,2)=1
      SHIFTXY(4,1)=1
      SHIFTXY(4,2)=1
      WRITE(WINDBG%FHNDL,*) 'LSAVE_INTERP_ARRAY=', LSAVE_INTERP_ARRAY
      IF (LSAVE_INTERP_ARRAY) THEN
        CALL LOAD_INTERP_ARRAY(FileSave, success)
        WRITE(WINDBG%FHNDL,*) 'success=', success
        IF (success .eqv. .TRUE.) RETURN
      END IF
      CF_IX=0
      CF_IY=0
      CF_COEFF=0
      WRITE(WINDBG%FHNDL,*) 'min(lon)=', minval(TheInfo % LON)
      WRITE(WINDBG%FHNDL,*) 'max(lon)=', maxval(TheInfo % LON)
      WRITE(WINDBG%FHNDL,*) 'min(lat)=', minval(TheInfo % LAT)
      WRITE(WINDBG%FHNDL,*) 'max(lat)=', maxval(TheInfo % LAT)
      WRITE(WINDBG%FHNDL,*) 'MNP_WIND=', MNP_WIND
      WRITE(WINDBG%FHNDL,*) 'START BIG SLOW LOOP'
      FLUSH(WINDBG%FHNDL)
      DO I = 1, MNP_WIND
        IF (I .eq. 1) THEN
          IXs=1
          IYs=1
        ELSE
          IXs=CF_IX(I-1)
          IYs=CF_IX(I-1)
        END IF
        eX=XP_WIND(I)
        eY=YP_WIND(I)
        CALL COMPUTE_SINGLE_INTERPOLATION_INFO(TheInfo, EXTRAPOLATION_ALLOWED_WIND, eX, eY, eCF_IX, eCF_IY, eCF_COEFF, EXTRAPO_OUT)
        CF_IX(I) = eCF_IX
        CF_IY(I) = eCF_IY
        CF_COEFF(:,I) = eCF_COEFF
        IF (EXTRAPO_OUT .eqv. .TRUE.) THEN
          nbExtrapolation = nbExtrapolation + 1
        END IF
	WRITE(WINDBG%FHNDL,*) 'I(', MNP_WIND,')=', I
	FLUSH(WINDBG%FHNDL)
      END DO
      IF (LSAVE_INTERP_ARRAY) THEN
        CALL SAVE_INTERP_ARRAY(FileSave)
      END IF
      IF (EXTRAPOLATION_ALLOWED_WIND .eqv. .TRUE.) THEN
        WRITE(WINDBG%FHNDL,*) ' nbExtrapolation=', nbExtrapolation
        WRITE(WINDBG%FHNDL,*) ' MaxMinDist=', sqrt(MaxMinDist)
      END IF
      WRITE(WINDBG%FHNDL,*) ' done interp calcs'
      FLUSH(WINDBG%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_CF_TIME_INDEX(eVAR, REC1, REC2, w1, w2)
      ! For given wwm_time and wind_time return records to get and weights for time
      ! interpolation F(wwm_time)=F(rec1)*w1 + F(rec2)*w2
      !
      USE DATAPOOL
      IMPLICIT NONE
      TYPE(VAR_NETCDF_CF), intent(in)           :: eVAR
      REAL(rkind), INTENT(OUT)            :: w1, w2
      INTEGER, INTENT(OUT)                :: REC1, REC2
      REAL(rkind) :: eTime1, eTime2
      INTEGER  :: iTime
 
      DO iTime=2,eVAR % nbTime
        eTime1=eVAR % ListTime(iTime-1)
        eTime2=eVAR % ListTime(iTime)
        IF ((eTime1 .le. MAIN%TMJD).and.(MAIN%TMJD .le. eTime2)) THEN
          REC2=iTime
          REC1=iTime-1
          w2=(MAIN % TMJD - eTime1)/(eTime2-eTime1)
          w1=(eTime2 - MAIN % TMJD)/(eTime2-eTime1)
          RETURN
        END IF
      END DO
      WRITE(WINDBG%FHNDL,*) 'Time error in wind for CF'
      WRITE(WINDBG%FHNDL,*) 'MAIN % TMJD=', MAIN%TMJD
      WRITE(WINDBG%FHNDL,*) 'min(wind_time_mjd)=', minval(eVAR % ListTime)
      WRITE(WINDBG%FHNDL,*) 'max(wind_time_mjd)=', maxval(eVAR % ListTime)
      FLUSH(WINDBG%FHNDL)
      CALL WWM_ABORT('Error in CF wind forcing time setup')
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SYNCHRONIZE_WIND_TIME_IFILE_IT
      USE DATAPOOL
      IMPLICIT NONE
      integer IPROC, eInt(1)
#ifdef MPI_PARALL_GRID
      IF (.NOT. MULTIPLE_IN_WIND) THEN
        IF (myrank .eq. 0) THEN
          eInt(1)=NDT_WIND_ALL_FILES
          DO IPROC=2,nproc
            CALL MPI_SEND(eInt,1,itype, iProc-1, 811, comm, ierr)
          END DO
          DO IPROC=2,nproc
            CALL MPI_SEND(WIND_TIME_ALL_FILES,NDT_WIND_ALL_FILES,rtype, iProc-1, 812, comm, ierr)
            CALL MPI_SEND(WIND_TIME_IFILE,NDT_WIND_ALL_FILES,itype, iProc-1, 813, comm, ierr)
            CALL MPI_SEND(WIND_TIME_IT,NDT_WIND_ALL_FILES,itype, iProc-1, 814, comm, ierr)
          END DO
        ELSE
          CALL MPI_RECV(eInt,1,itype, 0, 811, comm, istatus, ierr)
          NDT_WIND_ALL_FILES=eInt(1)
          ALLOCATE(WIND_TIME_ALL_FILES(NDT_WIND_ALL_FILES), WIND_TIME_IFILE(NDT_WIND_ALL_FILES), WIND_TIME_IT(NDT_WIND_ALL_FILES), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 3')
          CALL MPI_RECV(WIND_TIME_ALL_FILES,NDT_WIND_ALL_FILES,rtype, 0, 812, comm, istatus, ierr)
          CALL MPI_RECV(WIND_TIME_IFILE,NDT_WIND_ALL_FILES,rtype, 0, 813, comm, istatus, ierr)
          CALL MPI_RECV(WIND_TIME_IT,NDT_WIND_ALL_FILES,rtype, 0, 814, comm, istatus, ierr)
        END IF
      END IF
#endif
      SEWI%DELT = ( WIND_TIME_ALL_FILES(2) - WIND_TIME_ALL_FILES(1) ) * DAY2SEC
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef NCDF
      SUBROUTINE INIT_NETCDF_DWD
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
      INTEGER :: IT, IFILE
      INTEGER :: ILON_ID, ILAT_ID, ITIME_ID, I, J, COUNTER
      REAL(rkind)  :: DTMP
      REAL(rkind), ALLOCATABLE :: WIND_TIME(:)
      character ( len = 20 ) chrtmp
      character ( len = 15 ) chrdate

      integer, dimension(nf90_max_var_dims) :: dimIDs
      character (len = *), parameter :: CallFct="INIT_NETCDF_DWD"
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
# endif
        OPEN(WIN%FHNDL,FILE=WIN%FNAME,STATUS='OLD')
!
! count number of netcdf files in list ...
!
        NUM_NETCDF_FILES = 0
        DO
          READ( WIN%FHNDL, *, IOSTAT = ISTAT )
          IF ( ISTAT /= 0 ) EXIT
          NUM_NETCDF_FILES = NUM_NETCDF_FILES + 1
        END DO
        REWIND (WIN%FHNDL)

        ALLOCATE(NETCDF_FILE_NAMES(NUM_NETCDF_FILES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 3')

        DO IT = 1, NUM_NETCDF_FILES
          READ( WIN%FHNDL, *) NETCDF_FILE_NAMES(IT)
!          WRITE(WINDBG%FHNDL,*) IT, NETCDF_FILE_NAMES(IT)
        END DO
        CLOSE (WIN%FHNDL)
!
! check number of time steps in netcdf file ... it is assumed that all files have the same ammount of time steps ...
!
        CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(1)))
        ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(1), NF90_NOWRITE, WIND_NCID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)

        ISTAT = nf90_inq_varid(WIND_NCID, 'initial_time0_encoded', ITIME_ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, ITIME_ID, dimids = dimids)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = NDT_WIND_FILE)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, ISTAT)
!
! check dimensions in the netcdf ... again it is assumed that this is not changing for all files ...
!
        ISTAT = nf90_inq_varid(WIND_NCID, 'g0_lon_2', ILON_ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, ILON_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = NDX_WIND)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, ISTAT)

        ISTAT = nf90_inq_varid(WIND_NCID, 'g0_lat_1', ILAT_ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, ILAT_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = NDY_WIND)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, ISTAT)

        ALLOCATE (COORD_WIND_X(NDX_WIND), COORD_WIND_Y(NDY_WIND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 4')
!
! read coordinates from files ....
!
        ISTAT = NF90_GET_VAR(WIND_NCID, ILON_ID, COORD_WIND_X)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)

        ISTAT = NF90_GET_VAR(WIND_NCID, ILAT_ID, COORD_WIND_Y)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, ISTAT)
!
! estimate offset ...
!
        OFFSET_X_WIND = MINVAL(COORD_WIND_X)
        OFFSET_Y_WIND = MINVAL(COORD_WIND_Y)
!
! resolution ...
!
        DX_WIND  = ABS(MAXVAL(COORD_WIND_X)-MINVAL(COORD_WIND_X))/(NDX_WIND-1)
        DY_WIND  = ABS(MAXVAL(COORD_WIND_Y)-MINVAL(COORD_WIND_Y))/(NDY_WIND-1)
!
! close netcdf file ...
!
        ISTAT = NF90_CLOSE(WIND_NCID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, ISTAT)
!
! total number of time steps ... in all files
!
        NDT_WIND_ALL_FILES = NDT_WIND_FILE * NUM_NETCDF_FILES

        ALLOCATE (WIND_TIME(NDT_WIND_FILE), WIND_TIME_ALL_FILES(NDT_WIND_ALL_FILES), WIND_TIME_IFILE(NDT_WIND_ALL_FILES), WIND_TIME_IT(NDT_WIND_ALL_FILES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 5')
!
! read all time steps in the proper format and transform in wwm time line
!
        DO IFILE = 1, NUM_NETCDF_FILES
          CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(IFILE)))
          ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WIND_NCID)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, ISTAT)

          ISTAT = NF90_GET_VAR(WIND_NCID, ITIME_ID, WIND_TIME)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 15, ISTAT)

          DO IT = 1, NDT_WIND_FILE
            WRITE (3001, *) WIND_TIME(IT)
          END DO
          REWIND(3001)
          DO IT = 1, NDT_WIND_FILE
            READ (3001, '(A20)') CHRTMP
            ! YYYYMMDD
            CHRDATE(1:8) = CHRTMP(4:12)
            CHRDATE(9:9) = '.'
            ! HH
            CHRDATE(10:11) = CHRTMP(12:13)
            ! MMSS
            CHRDATE(12:15) = '0000' ! Construct propper format YYYYMMDDHHMMSS character .... len = 15 ... :)
            CALL CT2MJD(CHRDATE,DTMP) !
            WIND_TIME(IT) = DTMP    ! Double time with respect to 19000101.000000
          END DO
          CLOSE(3001)
          DO IT = 1, NDT_WIND_FILE
            WIND_TIME_ALL_FILES(IT+(IFILE-1)*NDT_WIND_FILE) = WIND_TIME(IT)
            WIND_TIME_IFILE(IT+(IFILE-1)*NDT_WIND_FILE) = IFILE
            WIND_TIME_IT   (IT+(IFILE-1)*NDT_WIND_FILE) = IT
          END DO
        END DO ! IFILE
        IF (LWRITE_ORIG_WIND) THEN
          WRITE (3010, '(I10)') 0
          WRITE (3010, '(I10)') NDX_WIND * NDY_WIND
          COUNTER = 0
          DO I = 1, NDY_WIND
            DO J = 1, NDX_WIND
              WRITE (3010, '(I10,3F15.4)') COUNTER, OFFSET_X_WIND+(J-1)*DX_WIND ,OFFSET_Y_WIND+(I-1)*DY_WIND , 0.0
              COUNTER = COUNTER + 1
            END DO
          END DO
          WRITE (3010, *) (NDX_WIND-1)*(NDY_WIND-1)*2
          DO J = 0, NDY_WIND-2
            DO I = 0, NDX_WIND-2
              WRITE (3010, '(5I10)')  I+J*NDX_WIND           , NDX_WIND+I+J* NDX_WIND, NDX_WIND+I+1+J*NDX_WIND, 0, 0
              WRITE (3010, '(5I10)')  NDX_WIND+I+1+J*NDX_WIND, I+1+J*NDX_WIND        , I+J*NDX_WIND           , 0, 0
            END DO
          END DO
          OPEN(3011, FILE  = 'ergwindorig.bin', FORM = 'UNFORMATTED')
        END IF
        ALLOCATE (WIND_X(NDX_WIND,NDY_WIND), WIND_Y(NDX_WIND,NDY_WIND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 9')
# ifdef MPI_PARALL_GRID
      END IF
# endif
      CALL SYNCHRONIZE_WIND_TIME_IFILE_IT
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_DWD(IFILE, IT, eField)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
!
!        READS WIND_Y, WIND_X and PRESSURE from a given NCID within one DWD file
!
      INTEGER, INTENT(IN) :: IFILE, IT
      REAL(rkind), intent(inout) :: eField(MNP,2)
      REAL(rkind)                :: Vtotal1(MNP_WIND)
      REAL(rkind)                :: Vtotal2(MNP_WIND)
      REAL(rkind)                :: Vlocal(MNP)
      character (len = *), parameter :: CallFct="READ_NETCDF_DWD"
      INTEGER             :: DWIND_X_ID, DWIND_Y_ID
      INTEGER             :: numLons, numLats, numTime, iy, counter, ip, i, j
      REAL(rkind),   ALLOCATABLE :: TMP(:,:)
      REAL(rkind), ALLOCATABLE   :: U(:), V(:), H(:)
      REAL(rkind), SAVE          :: TIME

      INTEGER, DIMENSION (nf90_max_var_dims) :: dimIDs
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
# endif
        IF (IFILE .GT. NUM_NETCDF_FILES) CALL WWM_ABORT('SOMETHING IS WRONG WE RUN OUT OF WIND TIME')

        CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(IFILE)))
        ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WIND_NCID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)

        ISTAT = nf90_inq_varid(WIND_NCID, 'V_GDS0_HTGL_13', DWIND_X_ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, DWIND_X_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = numLons)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(2), len = numLats)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(3), len = numTime)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, ISTAT)

        ISTAT = nf90_inq_varid(WIND_NCID, 'U_GDS0_HTGL_13', DWIND_Y_ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, DWIND_Y_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = numLons)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(2), len = numLats)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(3), len = numTime)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)

        ISTAT = NF90_GET_VAR(WIND_NCID, DWIND_X_ID, WIND_X,    start = (/ 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND, 1 /))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, ISTAT)

        ISTAT = NF90_GET_VAR(WIND_NCID, DWIND_Y_ID, WIND_Y,    start = (/ 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND, 1 /))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, ISTAT)

        IF (LINVERTY) THEN
          ALLOCATE(TMP(NDX_WIND,NDY_WIND), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 11')

          DO IY = 1, NDY_WIND
            TMP(:,NDY_WIND-(IY-1)) = wind_x(:,IY)
          END DO
          wind_x = TMP
          DO IY = 1, NDY_WIND
            TMP(:,NDY_WIND-(IY-1)) = wind_y(:,IY)
          END DO
          wind_y = TMP
          DO IY = 1, NDY_WIND
            TMP(:,NDY_WIND-(IY-1)) = atmo_press(:,IY)
          END DO
          atmo_press = TMP
          DEALLOCATE(TMP)
        END IF

        IF (LWRITE_ORIG_WIND) THEN
          ALLOCATE(U(NDX_WIND*NDY_WIND), V(NDX_WIND*NDY_WIND), H(NDX_WIND*NDY_WIND), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 15')
          COUNTER = 1
          DO J = 1, NDY_WIND
            DO I = 1, NDX_WIND
              U(COUNTER) = WIND_X(I,J)
              V(COUNTER) = WIND_Y(I,J)
              IF (ABS(U(COUNTER)) .GT. 1000.) U(COUNTER) = 0.
              IF (ABS(V(COUNTER)) .GT. 1000.) V(COUNTER) = 0.
              H(COUNTER) = SQRT((U(COUNTER)**2.+V(COUNTER)**2.))
              COUNTER = COUNTER + 1
            END DO
          END DO

          TIME = TIME + 1.
          WRITE(3011) TIME
          WRITE(3011) (U(IP), V(IP), H(IP), IP = 1, numLons*numLats)
          DEALLOCATE(U,V,H)
        END IF
        ISTAT = NF90_CLOSE(WIND_NCID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, ISTAT)

        CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_X,Vtotal1)
        CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_Y,Vtotal2)
# ifdef MPI_PARALL_GRID
      END IF
# endif
# ifdef MPI_PARALL_GRID
     IF (MULTIPLE_IN_WIND) THEN
       eField(:,1)=Vtotal1
       eField(:,2)=Vtotal2
     ELSE
       CALL SCATTER_ONED_ARRAY(Vtotal1, Vlocal)
       eField(:,1)=Vlocal
       CALL SCATTER_ONED_ARRAY(Vtotal2, Vlocal)
       eField(:,2)=Vlocal
     END IF
# else
     eField(:,1)=Vtotal1
     eField(:,2)=Vtotal2
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_NETCDF_CRFS
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
!for data description consult ftp://nomads.ncdc.noaa.gov/CFSR/HP_time_series/200307/wnd10m.l.gdas.200307.grb2.inv
      INTEGER :: IT, IFILE
      INTEGER :: ILON_ID, ILAT_ID, ITIME_ID, I, J, COUNTER
      REAL(rkind)   :: START_TIME
      REAL(rkind) , ALLOCATABLE :: WIND_TIME(:), WIND_TIME_NETCDF(:)
      character ( len = 15 ) chrdate
      character ( len = 40 ) beginn_time
      CHARACTER(LEN=15) :: eTimeStr
      character (len = *), parameter :: CallFct="INIT_NETCDF_CRFS"
      integer, dimension(nf90_max_var_dims) :: dimIDs
      logical PREF_ANALYZED
      integer idx
      REAL(rkind) :: ePresTime, eNewTime
      WRITE(WINDBG%FHNDL,*) 'Begin INIT_NETCDF_CRFS'

# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
# endif
        OPEN(WIN%FHNDL,FILE=WIN%FNAME,STATUS='OLD',IOSTAT = ISTAT)
!
! count number of netcdf files in list ...
! CRFS has analyzed at time 0 and forecast at times +1, +2, +3, +4, +5 +6
! Ordering in the file is
! --analyzed 0
! --fcst +1
! --fcst +2
! --fcst +3
! --fcst +4
! --fcst +5
! --fcst +6
! So, it goes by blocks of 7.
! PREF_ANALYZED = TRUE.  For prefering analyzed fields when available
!                 FALSE  For using the analyzed field only at first step.
        PREF_ANALYZED=.FALSE.
        NUM_NETCDF_FILES = 0
        DO
          READ( WIN%FHNDL, *, IOSTAT = ISTAT )
          IF ( ISTAT /= 0 ) EXIT
          NUM_NETCDF_FILES = NUM_NETCDF_FILES + 1
        END DO
        REWIND (WIN%FHNDL)

        ALLOCATE(NETCDF_FILE_NAMES(NUM_NETCDF_FILES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 18')

        DO IT = 1, NUM_NETCDF_FILES
          READ( WIN%FHNDL, *) NETCDF_FILE_NAMES(IT)
          WRITE(WINDBG%FHNDL,*) IT, NETCDF_FILE_NAMES(IT)
        END DO
        CLOSE (WIN%FHNDL)
!
! check number of time steps in netcdf file ... it is assumed that all files have the same ammount of time steps ...
!
        CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(1)))
        ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(1), NF90_NOWRITE, WIND_NCID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)

        ISTAT = nf90_inq_varid(WIND_NCID, 'time', ITIME_ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, ITIME_ID, dimids = dimids)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = NDT_WIND_FILE)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, ISTAT)

        ISTAT = nf90_get_att(WIND_NCID, ITIME_ID, 'units', beginn_time)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)

        CHRDATE(1:4) = beginn_time(13:17)
        CHRDATE(5:6) = beginn_time(18:19)
        CHRDATE(7:8) = beginn_time(21:22)
        CHRDATE(9:9)   = '.'
        CHRDATE(10:15)= '000000'
        CALL CT2MJD(CHRDATE,START_TIME)
!
! check dimensions in the netcdf ... again it is assumed that this is not changing for all files ...
!
        ISTAT = nf90_inq_varid(WIND_NCID, 'lon', ILON_ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, ILON_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = NDX_WIND)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, ISTAT)

        ISTAT = nf90_inq_varid(WIND_NCID, 'lat', ILAT_ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, ILAT_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = NDY_WIND)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)

        ALLOCATE (DCOORD_WIND_X(NDX_WIND), DCOORD_WIND_Y(NDY_WIND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 19')
!
! read cooridantes from files ....
!
        ISTAT = NF90_GET_VAR(WIND_NCID, ILON_ID, DCOORD_WIND_X)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, ISTAT)

        ISTAT = NF90_GET_VAR(WIND_NCID, ILAT_ID, DCOORD_WIND_Y)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, ISTAT)

        DCOORD_WIND_Y = DCOORD_WIND_Y !+ 90.0_rkind

        DO I = 1, NDX_WIND
          DCOORD_WIND_X(I) = DCOORD_WIND_X(I) - 180.0_rkind
        END DO
!
! estimate offset ...
!
        OFFSET_X_WIND = MINVAL(DCOORD_WIND_X)
        OFFSET_Y_WIND = MINVAL(DCOORD_WIND_Y)
!
! resolution ...
!
        DX_WIND  = ABS(MAXVAL(DCOORD_WIND_X)-MINVAL(DCOORD_WIND_X))/(NDX_WIND-1)
        DY_WIND  = ABS(MAXVAL(DCOORD_WIND_Y)-MINVAL(DCOORD_WIND_Y))/(NDY_WIND-1)
!
! close netcdf file ...
!
        ISTAT = NF90_CLOSE(WIND_NCID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, ISTAT)
!
! total number of time steps ... in all files
!
        ALLOCATE (WIND_X(NDX_WIND,NDY_WIND), WIND_Y(NDX_WIND,NDY_WIND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 35')

        ALLOCATE (WIND_TIME(NDT_WIND_FILE),WIND_TIME_NETCDF(NDT_WIND_FILE), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 20')
        WRITE(WINDBG%FHNDL,*) 'NDT_WIND_FILE=', NDT_WIND_FILE
        WRITE(WINDBG%FHNDL,*) 'START_TIME=', START_TIME
!
! read all time steps in the proper format and transform in wwm time line
!
        NDT_WIND_ALL_FILES=0
        ePresTime=-100000000
        DO IFILE = 1, NUM_NETCDF_FILES
          CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(IFILE)))
          ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WIND_NCID)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 15, ISTAT)

          ISTAT = NF90_GET_VAR(WIND_NCID, ITIME_ID, WIND_TIME_NETCDF)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 16, ISTAT)

          WRITE(WINDBG%FHNDL,*) 'IFILE=', IFILE
          DO IT = 1, NDT_WIND_FILE
            eNewTime=START_TIME+WIND_TIME_NETCDF(IT)*3600. * SEC2DAY
            IF (eNewTime .gt. ePresTime + THR8) THEN
              NDT_WIND_ALL_FILES=NDT_WIND_ALL_FILES + 1
              ePresTime=eNewTime
            END IF
          END DO
        END DO
        WRITE(WINDBG%FHNDL,*) 'NDT_WIND_ALL_FILES=', NDT_WIND_ALL_FILES
        ALLOCATE (WIND_TIME_ALL_FILES(NDT_WIND_ALL_FILES), WIND_TIME_IFILE(NDT_WIND_ALL_FILES), WIND_TIME_IT(NDT_WIND_ALL_FILES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 21')

        idx=0
        ePresTime=-100000000
        DO IFILE = 1, NUM_NETCDF_FILES
          CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(IFILE)))
          ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WIND_NCID)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 17, ISTAT)

          ISTAT = NF90_GET_VAR(WIND_NCID, ITIME_ID, WIND_TIME_NETCDF)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 18, ISTAT)

          WRITE(WINDBG%FHNDL,*) 'IFILE=', IFILE
          DO IT = 1, NDT_WIND_FILE
            eNewTime=START_TIME+WIND_TIME_NETCDF(IT)*3600. * SEC2DAY
            CALL MJD2CT(eNewTime,eTimeStr)
            IF (PREF_ANALYZED) THEN
              IF (eNewTime .gt. ePresTime + THR8) THEN
                ePresTime=eNewTime
                idx=idx+1
                WIND_TIME_ALL_FILES(idx) = eNewTime
                WIND_TIME_IFILE(idx) = IFILE
                WIND_TIME_IT(idx) = IT
                WRITE(WINDBG%FHNDL,110) IT, idx, WIND_TIME_NETCDF(IT) * 3600. * SEC2DAY, eNewTime, eTimeStr
              ENDIF
            ELSE
              IF (eNewTime .gt. ePresTime + THR8) THEN
                ePresTime=eNewTime
                idx=idx+1
              ENDIF
              WIND_TIME_ALL_FILES(idx) = eNewTime
              WIND_TIME_IFILE(idx) = IFILE
              WIND_TIME_IT(idx) = IT
              WRITE(WINDBG%FHNDL,110) IT, idx, WIND_TIME_NETCDF(IT) * 3600. * SEC2DAY, eNewTime, eTimeStr
            ENDIF
          END DO
        END DO ! IFILE
110     FORMAT (I4, ' ', I4, ' ', F15.3, ' ', F15.3, ' ', a15)

        IF (LWRITE_ORIG_WIND) THEN
          WRITE (3010, '(I10)') 0
          WRITE (3010, '(I10)') NDX_WIND * NDY_WIND
          COUNTER = 0
          DO I = 1, NDY_WIND
            DO J = 1, NDX_WIND
              WRITE (3010, '(I10,3F15.4)') COUNTER, OFFSET_X_WIND+(J-1)*DX_WIND ,OFFSET_Y_WIND+(I-1)*DY_WIND , 0.0
              COUNTER = COUNTER + 1
            END DO
          END DO
          WRITE (3010, *) (NDX_WIND-1)*(NDY_WIND-1)*2
          DO J = 0, NDY_WIND-2
            DO I = 0, NDX_WIND-2
              WRITE (3010, '(5I10)')  I+J*NDX_WIND           , NDX_WIND+I+J* NDX_WIND, NDX_WIND+I+1+J*NDX_WIND, 0, 0
              WRITE (3010, '(5I10)')  NDX_WIND+I+1+J*NDX_WIND, I+1+J*NDX_WIND        , I+J*NDX_WIND           , 0, 0
            END DO
          END DO
          OPEN(3011, FILE  = 'ergwindorig.bin', FORM = 'UNFORMATTED')
        END IF
# ifdef MPI_PARALL_GRID
      END IF
# endif
      CALL SYNCHRONIZE_WIND_TIME_IFILE_IT
      WRITE(WINDBG%FHNDL,*) 'End INIT_NETCDF_CRFS'
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_NETCDF_NARR
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
!
!2do update ...
!for data description consult ftp://nomads.ncdc.noaa.gov/CFSR/HP_time_series/200307/wnd10m.l.gdas.200307.grb2.inv
!
      INTEGER :: IT, IFILE, II, IP
      INTEGER :: ILON_ID, ILAT_ID, ITIME_ID, I, J, COUNTER
      REAL(rkind)  :: START_TIME, OFFSET_TIME
      REAL(rkind), ALLOCATABLE :: WIND_TIME_NETCDF(:)
      character ( len = 4 ) ch4
      character ( len = 15 ) chrdate
      character (len = *), parameter :: CallFct="INIT_NETCDF_NARR"
      integer, dimension(nf90_max_var_dims) :: dimIDs
      integer, allocatable :: COUNTERMAT(:,:)
      integer, allocatable :: IMAT(:), JMAT(:)
      integer :: NbPoint, nbFail
      real(rkind) :: Wi(3), XPW(3), YPW(3)
      INTEGER NI(3)
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
# endif
        WRITE(WINDBG%FHNDL,*) 'Begin INIT_NETCDF_NARR'
        WRITE(WINDBG%FHNDL,*) 'MULTIPLE_IN_WIND=', MULTIPLE_IN_WIND
!
! I make the assumption that the year when the dataset beginns at the year indicated in the bouc section
!
        CHRDATE(1:4) = SEBO%BEGT(1:4)
        CHRDATE(5:6) = '01'
        CHRDATE(7:8) = '01'
        CHRDATE(9:9)   = '.'
        CHRDATE(10:15)= '000000'
        CALL CT2MJD(CHRDATE,START_TIME)

        OPEN(WIN%FHNDL,FILE=WIN%FNAME,STATUS='OLD',IOSTAT = ISTAT)
!
! count number of netcdf files in list ...
!
        NUM_NETCDF_FILES = 0
        DO
          READ( WIN%FHNDL, *, IOSTAT = ISTAT )
          IF ( ISTAT /= 0 ) EXIT
          NUM_NETCDF_FILES = NUM_NETCDF_FILES + 1
        END DO
        REWIND (WIN%FHNDL)
        ALLOCATE(NETCDF_FILE_NAMES(NUM_NETCDF_FILES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 22')
        WRITE(WINDBG%FHNDL,*) 'NUM_NETCDF_FILES=', NUM_NETCDF_FILES

        DO IT = 1, NUM_NETCDF_FILES
          READ( WIN%FHNDL, *) NETCDF_FILE_NAMES(IT)
          WRITE(WINDBG%FHNDL,*) 'IT=', IT, 'file=', NETCDF_FILE_NAMES(IT)
        END DO
        CLOSE (WIN%FHNDL)
!
! check number of time steps in netcdf file ... it is assumed that all files have the same ammount of time steps ...
!
        CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(1)))
        ISTAT = NF90_OPEN(TRIM(NETCDF_FILE_NAMES(1)), NF90_NOWRITE, WINDX_NCID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)

        ISTAT = nf90_inq_varid(WINDX_NCID, 'time', ITIME_ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WINDX_NCID, ITIME_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, ISTAT)

        ISTAT = nf90_inquire_dimension(WINDX_NCID, dimIDs(1), len = NDT_WIND_FILE)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, ISTAT)

        WRITE(WINDBG%FHNDL,*) 'NDT_WIND_FILE=', NDT_WIND_FILE

        CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(2)))
        ISTAT = NF90_OPEN(TRIM(NETCDF_FILE_NAMES(2)), NF90_NOWRITE, WINDY_NCID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)
!
! check dimensions in the netcdf ... again it is assumed that this is not changing for all files ...
!
        ISTAT = nf90_inq_varid(WINDX_NCID, 'lon', ILON_ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WINDX_NCID, ILON_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, ISTAT)

        ISTAT = nf90_inquire_dimension(WINDX_NCID, dimIDs(1), len = NDX_WIND)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, ISTAT)

        ISTAT = nf90_inquire_dimension(WINDX_NCID, dimIDs(2), len = NDY_WIND)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, ISTAT)

        ISTAT = nf90_inq_varid(WINDY_NCID, 'lat', ILAT_ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WINDY_NCID, ILAT_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)

        ISTAT = nf90_inquire_dimension(WINDY_NCID, dimIDs(1), len = NDX_WIND)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, ISTAT)

        ISTAT = nf90_inquire_dimension(WINDY_NCID, dimIDs(2), len = NDY_WIND)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, ISTAT)

        ALLOCATE (DCOORD_WIND_X2(NDX_WIND,NDY_WIND), DCOORD_WIND_Y2(NDX_WIND,NDY_WIND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 23')
!
! read coordinates from files ....
!
        ISTAT = NF90_GET_VAR(WINDX_NCID, ILON_ID, DCOORD_WIND_X2)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, ISTAT)

        ISTAT = NF90_GET_VAR(WINDY_NCID, ILAT_ID, DCOORD_WIND_Y2)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 15, ISTAT)
!
! estimate offset ...
!
        OFFSET_X_WIND = MINVAL(DCOORD_WIND_X2)
        OFFSET_Y_WIND = MINVAL(DCOORD_WIND_Y2)
        WRITE(WINDBG%FHNDL,*) 'OFFSET_X_WIND=', OFFSET_X_WIND
        WRITE(WINDBG%FHNDL,*) 'OFFSET_Y_WIND=', OFFSET_Y_WIND
!
! close netcdf file ...
!
        ISTAT = NF90_CLOSE(WINDX_NCID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 16, ISTAT)

        ISTAT = NF90_CLOSE(WINDY_NCID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 17, ISTAT)
!
! total number of time steps ... in all files
!
        NDT_WIND_ALL_FILES = NDT_WIND_FILE * NUM_NETCDF_FILES/2
        WRITE(WINDBG%FHNDL,*) 'NDT_WIND_ALL_FILES=', NDT_WIND_ALL_FILES

        ALLOCATE (WIND_TIME_NETCDF(NDT_WIND_FILE), WIND_TIME_ALL_FILES(NDT_WIND_ALL_FILES), WIND_TIME_IFILE(NDT_WIND_ALL_FILES), WIND_TIME_IT(NDT_WIND_ALL_FILES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 24')
!
! read all time steps in the proper format and transform in wwm time line
!
        DO IFILE = 1, NUM_NETCDF_FILES, 2
          CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(IFILE)))
          ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WINDX_NCID)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 18, ISTAT)

          ISTAT = NF90_GET_VAR(WINDX_NCID, ITIME_ID, WIND_TIME_NETCDF)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 19, ISTAT)

          ISTAT = NF90_CLOSE(WINDX_NCID)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 20, ISTAT)

          WIND_TIME_NETCDF = WIND_TIME_NETCDF/24 ! Transform to days ...
          OFFSET_TIME = MINVAL(WIND_TIME_NETCDF) - START_TIME  ! in this dataset the time start at 1800 in WWM it start at 1900
          WIND_TIME_NETCDF = WIND_TIME_NETCDF - OFFSET_TIME ! Now in WWM timeline ...
          DO IT = 1, NDT_WIND_FILE
            WIND_TIME_ALL_FILES(IT+(IFILE-1)*NDT_WIND_FILE) = WIND_TIME_NETCDF(IT)
            WIND_TIME_IFILE(IT+(IFILE-1)*NDT_WIND_FILE) = IFILE
            WIND_TIME_IT(IT+(IFILE-1)*NDT_WIND_FILE) = IT
          END DO
          CH4 = CHRDATE(1:4)
          READ (CH4,'(i4)') I
          I = I + 1
          WRITE (11111,'(i4)') I
          REWIND(11111)
          READ(11111,*) ch4
          CHRDATE(1:4) = ch4
          CHRDATE(5:6) = '01'
          CHRDATE(7:8) = '01'
          CHRDATE(9:9)   = '.'
          CHRDATE(10:15)= '000000'
          CALL CT2MJD(CHRDATE,START_TIME)
          !WRITE(WINDBG%FHNDL,*) CHRDATE, START_TIME
        END DO ! IFILE
        !
        ! Now the geographic interpolation
        !
        NE_WIND = (NDX_WIND-1)*(NDY_WIND-1)*2
        NP_WIND =  NDX_WIND*NDY_WIND

        IF (LWRITE_ORIG_WIND) THEN
          WRITE (3010, '(I10)') 0
          WRITE (3010, '(I10)') NP_WIND
        END IF
        COUNTER = 0
        NbPoint=NDX_WIND*NDY_WIND
        ALLOCATE(XYPWIND(2,NbPoint), IMAT(NbPoint), JMAT(NbPoint), COUNTERMAT(NDX_WIND,NDY_WIND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 28')
        DO I = 1, NDX_WIND
          DO J = 1, NDY_WIND
            IF (DCOORD_WIND_X2(I,J) .GT. 0.) THEN
! Transformed to a unified domain extending below -180.
              DCOORD_WIND_X2(I,J) = -1 * DCOORD_WIND_X2(I,J) - (180.-DCOORD_WIND_X2(I,J)) * 2
            END IF
            COUNTER = COUNTER + 1
            XYPWIND(1,COUNTER) = DCOORD_WIND_X2(I,J)
            XYPWIND(2,COUNTER) = DCOORD_WIND_Y2(I,J)
            IMAT(COUNTER)=I
            JMAT(COUNTER)=J
            COUNTERMAT(I,J)=COUNTER
            IF (LWRITE_ORIG_WIND) WRITE (3010, '(I10,3F15.4)') COUNTER, XYPWIND(1,COUNTER), XYPWIND(2,COUNTER), 0.0
          END DO
        END DO

        IF (LWRITE_ORIG_WIND) WRITE (3010, *) NE_WIND
        NE_WIND = (NDX_WIND-1)*(NDY_WIND-1)*2
        ALLOCATE(INE_WIND(3,NE_WIND), UWND_NARR(NDX_WIND*NDY_WIND), VWND_NARR(NDX_WIND*NDY_WIND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 30')

        INE_WIND = 0
        UWND_NARR = 0.
        VWND_NARR = 0.

        II = 0
        DO I = 1, NDX_WIND-1
          DO J = 1, NDY_WIND-1
            II = II + 1
            INE_WIND(1,II) = COUNTERMAT(I,J)
            INE_WIND(2,II) = COUNTERMAT(I+1,J)
            INE_WIND(3,II) = COUNTERMAT(I,J+1)
            II = II + 1
            INE_WIND(1,II) = COUNTERMAT(I+1,J+1)
            INE_WIND(2,II) = COUNTERMAT(I,J+1)
            INE_WIND(3,II) = COUNTERMAT(I+1,J)
          END DO
        END DO

        IF (LWRITE_ORIG_WIND) OPEN(3011, FILE  = 'ergwindorig.bin', FORM = 'UNFORMATTED')
        ALLOCATE(WIND_ELE(MNP_WIND), WI_NARR(MNP_WIND, 3), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 32')
        WIND_ELE = 0
        WI_NARR = 0
        nbFail=0
        DO IP = 1, MNP_WIND
          CALL FIND_ELE_WIND( NE_WIND, NP_WIND, INE_WIND, XYPWIND, XP_WIND(IP), YP_WIND(IP), WIND_ELE(IP))
          IF (WIND_ELE(IP) .eq. 0) THEN
            WRITE(WINDBG%FHNDL,*) 'POINT OF THE MESH IS OUT OF THE WIND FIELD', IP, XP(IP), YP(IP)
            nbFail=nbFail+1
          ELSE
            NI=INE_WIND(:,WIND_ELE(IP))
            XPW=XYPWIND(1,NI)
            YPW=XYPWIND(2,NI)
            CALL INTELEMENT_COEF(XPW, YPW,XP_WIND(IP),YP_WIND(IP),Wi)
            WI_NARR(IP,:)=Wi
!           WRITE(WINDBG%FHNDL,*) 'IP=', MNP, ' sumWi=', sum(Wi)
!           WRITE(WINDBG%FHNDL,*) 'IP=', MNP, ' minW=', minval(Wi), ' maxW=', maxval(Wi)
          ENDIF 
        END DO
        WRITE(WINDBG%FHNDL,*) 'MNP_WIND=', MNP_WIND, ' nbFail=', nbFail
        WRITE(WINDBG%FHNDL,*) 'NDX_WIND=', NDX_WIND
        WRITE(WINDBG%FHNDL,*) 'NDY_WIND=', NDY_WIND
        DEALLOCATE(IMAT, JMAT, COUNTERMAT)
        ALLOCATE (WIND_X4(NDX_WIND,NDY_WIND), WIND_Y4(NDX_WIND,NDY_WIND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 41')
        CALL SYNCHRONIZE_WIND_TIME_IFILE_IT
# ifdef MPI_PARALL_GRID
      END IF
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_CRFS(IFILE, IT, eField)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
!
!     READS WIND_Y, WIND_X and PRESSURE from a given NCID within one CRFS file
!
      INTEGER, INTENT(IN) :: IFILE, IT
      REAL(rkind), intent(inout) :: eField(MNP,2)
      character (len = *), parameter :: CallFct="READ_NETCDF_CRFS"

      INTEGER             :: DWIND_X_ID, DWIND_Y_ID
      INTEGER             :: numLons, numLats, numTime, numHeights, iy, counter, ip, i, j, ix
      REAL(rkind),   ALLOCATABLE :: TMP(:,:)
      REAL(rkind), ALLOCATABLE   :: U(:), V(:)
      REAL(rkind), SAVE          :: TIME
      REAL(rkind)                :: Vtotal1(MNP_WIND)
      REAL(rkind)                :: Vtotal2(MNP_WIND)
      REAL(rkind)                :: Vlocal(MNP)
      INTEGER, DIMENSION (nf90_max_var_dims) :: dimIDs
      WRITE(WINDBG%FHNDL,*) 'READ_NETCDF_CRFS IFILE=', IFILE, ' IT=', IT
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
# endif
        IF (IFILE .GT. NUM_NETCDF_FILES) CALL WWM_ABORT('SOMETHING IS WRONG WE RUN OUT OF WIND TIME')

        CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(IFILE)))
        ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WIND_NCID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)

        ISTAT = nf90_inq_varid(WIND_NCID, '10u', DWIND_X_ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, DWIND_X_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = numLons)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(2), len = numLats)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(3), len = numHeights)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(4), len = numTime)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, ISTAT)

        ISTAT = nf90_inq_varid(WIND_NCID, '10v', DWIND_Y_ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, DWIND_Y_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = numLons)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(2), len = numLats)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(3), len = numHeights)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(4), len = numTime)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, ISTAT)


        ISTAT = NF90_GET_VAR(WIND_NCID, DWIND_X_ID, WIND_X, start = (/ 1, 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND,1,1 /))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, ISTAT)

        ISTAT = NF90_GET_VAR(WIND_NCID, DWIND_Y_ID, WIND_Y, start = (/ 1, 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND,1,1 /))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 15, ISTAT)

        ISTAT = NF90_CLOSE(WIND_NCID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 16, ISTAT)

        CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_X,Vtotal1)

        CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_Y,Vtotal2)
# ifdef MPI_PARALL_GRID
      END IF
# endif
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND) THEN
        eField(:,1)=Vtotal1
        eField(:,2)=Vtotal2
      ELSE
        CALL SCATTER_ONED_ARRAY(Vtotal1, Vlocal)
        eField(:,1)=Vlocal
        CALL SCATTER_ONED_ARRAY(Vtotal2, Vlocal)
        eField(:,2)=Vlocal
      END IF
# else
      eField(:,1)=Vtotal1
      eField(:,2)=Vtotal2
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_NARR(IFILE, IT, eField)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
!
!     READS WIND_Y, WIND_X and PRESSURE from a given NCID within one NARR file
!
      INTEGER, INTENT(IN) :: IFILE, IT
      REAL(rkind), intent(out) :: eField(MNP,2)

      INTEGER             :: DWIND_X_ID, DWIND_Y_ID
      INTEGER             :: numLons, numLats, counter, ip, i, j, ix
      character (len = *), parameter :: CallFct="READ_NETCDF_NARR"

      REAL(rkind),   ALLOCATABLE :: TMP(:,:)
      REAL(rkind),SAVE           :: TIME
      REAL(rkind)                :: scale_factor
      REAL(rkind) :: sumWi, eF1, eF2
      REAL(rkind) :: Vtotal1(MNP_WIND)
      REAL(rkind) :: Vtotal2(MNP_WIND)
      REAL(rkind) :: Vlocal(MNP)
      INTEGER, DIMENSION (nf90_max_var_dims) :: dimIDs
      real(rkind) :: ErrorCoord, XPinterp, YPinterp
      INTEGER IEwind
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
# endif
        IF (2*IFILE .GT. NUM_NETCDF_FILES) THEN
          CALL WWM_ABORT('NARR ERROR: Not enough files')
        END IF
        CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(2*IFILE-1)))
        ISTAT = NF90_OPEN(TRIM(NETCDF_FILE_NAMES(2*IFILE-1)), NF90_NOWRITE, WINDX_NCID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)

        ISTAT = nf90_inq_varid(WINDX_NCID, 'uwnd', DWIND_X_ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WINDX_NCID, DWIND_X_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, ISTAT)

        ISTAT = nf90_inquire_dimension(WINDX_NCID, dimIDs(1), len = numLons)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, ISTAT)

        ISTAT = nf90_inquire_dimension(WINDX_NCID, dimIDs(2), len = numLats)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)

        ISTAT = nf90_get_att(WINDX_NCID, DWIND_X_ID, 'scale_factor', scale_factor)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, ISTAT)

        CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(2*IFILE)))
        ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(2*IFILE), NF90_NOWRITE, WINDY_NCID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, ISTAT)

        ISTAT = nf90_inq_varid(WINDY_NCID, 'vwnd', DWIND_Y_ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WINDY_NCID, DWIND_Y_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, ISTAT)

        ISTAT = nf90_inquire_dimension(WINDY_NCID, dimIDs(1), len = numLons)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)

        ISTAT = nf90_inquire_dimension(WINDY_NCID, dimIDs(2), len = numLats)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, ISTAT)

        ISTAT = NF90_GET_VAR(WINDX_NCID, DWIND_X_ID, WIND_X4, start = (/ 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND, 1 /))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, ISTAT)

        ISTAT = NF90_GET_VAR(WINDY_NCID, DWIND_Y_ID, WIND_Y4, start = (/ 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND, 1 /))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 15, ISTAT)

        ISTAT = NF90_CLOSE(WINDX_NCID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 16, ISTAT)

        ISTAT = NF90_CLOSE(WINDY_NCID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 17, ISTAT)

        IF (.FALSE.) THEN
          ALLOCATE(TMP(NDX_WIND,NDY_WIND), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 43')
          DO IX = 1, NDX_WIND
            tmp(NDX_WIND-(IX-1),:) = wind_x4(IX,:)
          END DO
          wind_x4 = tmp
          DO IX = 1, NDX_WIND
            tmp(NDX_WIND-(IX-1),:) = wind_y4(IX,:)
          END DO
          wind_y4 = tmp
!         DO IY = 1, NDY_WIND
!           tmp(:,NDY_WIND-(IY-1)) = atmo_press(:,IY)
!         END DO
!         atmo_press = tmp
          DEALLOCATE(TMP)
        END IF

        IF (.FALSE.) THEN
          ALLOCATE(TMP(NDX_WIND,NDY_WIND), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 44')
          DO IX = 1, NDX_WIND
            IF (IX .GT. NDX_WIND/2) THEN
              tmp(IX-NDX_WIND/2,:) = wind_x(IX,:)
            ELSE
              tmp(IX+NDX_WIND/2,:) = wind_x(IX,:)
            END IF
          END DO
          wind_x = tmp
          DO IX = 1, NDX_WIND
            IF (IX .GT. NDX_WIND/2) THEN
              tmp(IX-NDX_WIND/2,:) = wind_y(IX,:)
            ELSE
              tmp(IX+NDX_WIND/2,:) = wind_y(IX,:)
            END IF
          END DO
          wind_y = tmp
          DEALLOCATE(TMP)
        END IF

        COUNTER = 0
        DO I = 1, NDX_WIND
          DO J = 1, NDY_WIND
            COUNTER = COUNTER + 1
            IF (WIND_X4(I,J) .GT. 32766 .OR. WIND_Y4(I,J) .GT. 32766) THEN
              UWND_NARR(COUNTER) = 0.
              VWND_NARR(COUNTER) = 0.
            ELSE
              UWND_NARR(COUNTER) = WIND_X4(I,J) * scale_factor
              VWND_NARR(COUNTER) = WIND_Y4(I,J) * scale_factor
            END IF
          END DO
        END DO

        IF (LWRITE_ORIG_WIND) THEN
          TIME = TIME + 1.
          WRITE(3011) TIME
          WRITE(3011) (UWND_NARR(IP), VWND_NARR(IP), SQRT(UWND_NARR(IP)**2.+VWND_NARR(IP)**2.), IP = 1, NDX_WIND*NDY_WIND)
        END IF

        ErrorCoord=0
        DO IP = 1, MNP_WIND
          IEwind=WIND_ELE(IP)
          IF (IEwind .gt. 0) then 
            XPinterp=0
            YPinterp=0
            eF1=0
            eF2=0
            sumWi=0
            DO I=1,3
              XPinterp=XPinterp + WI_NARR(IP,I)*XYPWIND(1,INE_WIND(I,IEwind))
              YPinterp=YPinterp + WI_NARR(IP,I)*XYPWIND(2,INE_WIND(I,IEwind))
              eF1=eF1+WI_NARR(IP,I)*UWND_NARR(INE_WIND(I,IEwind))
              eF2=eF2+WI_NARR(IP,I)*VWND_NARR(INE_WIND(I,IEwind))
              sumWi=sumWi + WI_NARR(IP,I)
            END DO
            Vtotal1(IP)=eF1
            Vtotal2(IP)=eF2
            ErrorCoord=ErrorCoord + abs(XPinterp - XP_WIND(IP)) + abs(YPinterp - YP_WIND(IP)) + abs(sumWi - 1)
          ELSE
            Vtotal1(IP)=ZERO
            Vtotal2(IP)=ZERO
          ENDIF
        END DO
        WRITE(WINDBG%FHNDL,*) 'ErrorCoord=', ErrorCoord
# ifdef MPI_PARALL_GRID
      END IF
# endif
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND) THEN
        eField(:,1)=Vtotal1
        eField(:,2)=Vtotal2
      ELSE
        CALL SCATTER_ONED_ARRAY(Vtotal1, Vlocal)
        eField(:,1)=Vlocal
        CALL SCATTER_ONED_ARRAY(Vtotal2, Vlocal)
        eField(:,2)=Vlocal
      END IF
# else
      eField(:,1)=Vtotal1
      eField(:,2)=Vtotal2
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CHECK_WIND_TIME(nbtime_mjd, WIND_TIME_MJD)
      USE DATAPOOL, only : SEWI, WINDBG, rkind, THR, wwmerr
      IMPLICIT NONE
      integer, intent(in) :: nbtime_mjd
      real(rkind), intent(in) :: WIND_TIME_MJD(nbtime_mjd)
      CHARACTER(LEN=15) :: eTimeStr
      IF (SEWI%BMJD .LT. minval(WIND_TIME_MJD) - THR) THEN
        WRITE(WINDBG%FHNDL,*) 'END OF RUN'
        WRITE(WINDBG%FHNDL,*) 'WIND START TIME is outside CF wind_time range!'
        CALL MJD2CT(SEWI%BMJD,eTimeStr)
        WRITE(WINDBG%FHNDL,*) 'SEWI%BMJD=', SEWI%BMJD, ' date=', eTimeStr
        CALL MJD2CT(SEWI%EMJD,eTimeStr)
        WRITE(WINDBG%FHNDL,*) 'SEWI%EMJD=', SEWI%EMJD, ' date=', eTimeStr
        CALL MJD2CT(minval(WIND_TIME_MJD),eTimeStr)
        WRITE(WINDBG%FHNDL,*) 'min(WIND_TIME_MJD)=', minval(WIND_TIME_MJD), ' date=', eTimeStr
        CALL MJD2CT(maxval(WIND_TIME_MJD),eTimeStr)
        WRITE(WINDBG%FHNDL,*) 'max(WIND_TIME_MJD)=', maxval(WIND_TIME_MJD), ' date=', eTimeStr
        FLUSH(WINDBG%FHNDL)
        WRITE(wwmerr, *) 'Error in WIND_TIME_MJD 1, read ', TRIM(WINDBG%FNAME)
        CALL WWM_ABORT(wwmerr)
      END IF
      IF (SEWI%EMJD .GT. maxval(WIND_TIME_MJD) + THR) THEN
        WRITE(WINDBG%FHNDL,*) 'END OF RUN'
        WRITE(WINDBG%FHNDL,*) 'WIND END TIME is outside CF wind_time range!'
        CALL MJD2CT(SEWI%BMJD,eTimeStr)
        WRITE(WINDBG%FHNDL,*) 'SEWI%BMJD=', SEWI%BMJD, ' date=', eTimeStr
        CALL MJD2CT(SEWI%EMJD,eTimeStr)
        WRITE(WINDBG%FHNDL,*) 'SEWI%EMJD=', SEWI%EMJD, ' date=', eTimeStr
        CALL MJD2CT(minval(WIND_TIME_MJD),eTimeStr)
        WRITE(WINDBG%FHNDL,*) 'min(WIND_TIME_MJD)=', minval(WIND_TIME_MJD), ' date=', eTimeStr
        CALL MJD2CT(maxval(WIND_TIME_MJD),eTimeStr)
        WRITE(WINDBG%FHNDL,*) 'max(WIND_TIME_MJD)=', maxval(WIND_TIME_MJD), ' date=', eTimeStr
        FLUSH(WINDBG%FHNDL)
        WRITE(wwmerr, *) 'Error in WIND_TIME_MJD 2, read ', TRIM(WINDBG%FNAME)
        CALL WWM_ABORT(wwmerr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_INTERP_NETCDF_CF_WWM_WIND(RECORD_IN, outwind)
      USE NETCDF
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(in)                :: RECORD_IN
      REAL(rkind), INTENT(out)           :: outwind(MNP,2)
      REAL(rkind) :: varTotal(MNP_WIND,2), Vtotal(MNP_WIND), Vlocal(MNP)
      character (len = *), parameter :: CallFct="READ_INTERP_NETCDF_CF_WWM_WIND"
      INTEGER                            :: FID, ID
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
# endif
        CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(WIN%FNAME))
        ISTAT = NF90_OPEN(WIN%FNAME, NF90_NOWRITE, FID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)

        ISTAT = NF90_inq_varid(FID, 'Uwind', ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)

        ISTAT = NF90_GET_VAR(FID, ID, UWIND_FD, start = (/ 1, 1, RECORD_IN /), count = (/ NDX_WIND_FD, NDY_WIND_FD, 1 /))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, ISTAT)

        ISTAT = NF90_inq_varid(FID, 'Vwind', ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, ISTAT)

        ISTAT = NF90_GET_VAR(FID, ID, VWIND_FD, start = (/ 1, 1, RECORD_IN /), count = (/ NDX_WIND_FD, NDY_WIND_FD, 1 /))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)

        ISTAT = NF90_CLOSE(FID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, ISTAT)
        CALL KERNEL_INTERP_UV_WINDFD(varTotal)
# ifdef MPI_PARALL_GRID
      END IF
# endif
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND) THEN
        outwind=varTotal
      ELSE
        Vtotal=varTotal(:,1)
        CALL SCATTER_ONED_ARRAY(Vtotal, Vlocal)
        outwind(:,1)=Vlocal
        !
        Vtotal=varTotal(:,2)
        CALL SCATTER_ONED_ARRAY(Vtotal, Vlocal)
        outwind(:,2)=Vlocal
      END IF
# else
      outwind=varTotal
# endif
      END SUBROUTINE
!****************************************************************************
!*  CF_COMPLIANT WIND                                                       *
!*  This is the standard way to write netcdf data.                          *
!*  See                                                                     *
!* http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.6/cf-conventions.pdf *
!*  for details                                                             *
!****************************************************************************
      SUBROUTINE INIT_NETCDF_CF_WWM_WIND(eVAR)
      USE NETCDF
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER           :: fid, varid, dimids(2), dimidsB(3)
      integer nbChar
      TYPE(VAR_NETCDF_CF), intent(inout) :: eVAR
      REAL(rkind), ALLOCATABLE :: CF_LON(:,:), CF_LAT(:,:)
      character (len = *), parameter :: CallFct="INIT_NETCDF_CF"
      character (len=200) :: CoordString
      character (len=100) :: Xname, Yname, eStrUnitTime
      character (len=20) :: WindTimeStr
      real(rkind) :: ConvertToDay
      real(rkind) :: eTimeStart
      character(len=100) :: CHRERR
      integer posBlank, alen
      type(FD_FORCING_GRID) TheInfo
      integer IX, IY
      CALL INIT_DIRECT_NETCDF_CF(eVAR, MULTIPLE_IN_WIND, WIN%FNAME, "Uwind")
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
# endif
        ISTAT = nf90_open(WIN%FNAME, nf90_nowrite, fid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)

        ISTAT = nf90_inq_varid(fid, "Uwind", varid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)

        ISTAT = nf90_get_att(fid, varid, "coordinates", CoordString)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)
        alen=LEN_TRIM(CoordString)
        posBlank=INDEX(CoordString(1:alen), ' ')
        Xname=CoordString(1:posBlank-1)
        Yname=CoordString(posBlank+1:alen)
        WRITE(WINDBG%FHNDL,*) 'Xname=', TRIM(Xname)
        WRITE(WINDBG%FHNDL,*) 'Yname=', TRIM(Yname)
        FLUSH(WINDBG%FHNDL)

        ! Reading lontitude/latitude array

        ISTAT = nf90_inq_varid(fid, Xname, varid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, ISTAT)

        ISTAT = nf90_inquire_variable(fid, varid, dimids=dimids)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, ISTAT)

        ISTAT = nf90_inquire_dimension(fid, dimids(1), len=NDX_WIND_FD)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, ISTAT)

        ISTAT = nf90_inquire_dimension(fid, dimids(2), len=NDY_WIND_FD)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, ISTAT)

        WRITE(WINDBG%FHNDL,*) 'NDX_WIND_FD=', NDX_WIND_FD
        WRITE(WINDBG%FHNDL,*) 'NYX_WIND_FD=', NDY_WIND_FD
        FLUSH(WINDBG%FHNDL)

        allocate(CF_LON(NDX_WIND_FD, NDY_WIND_FD), CF_LAT(NDX_WIND_FD, NDY_WIND_FD), UWIND_FD(NDX_WIND_FD, NDY_WIND_FD), VWIND_FD(NDX_WIND_FD, NDY_WIND_FD), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 47')

        ISTAT = nf90_get_var(fid, varid, CF_LON)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, ISTAT)

        ISTAT = nf90_inq_varid(fid, Yname, varid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)

        ISTAT = nf90_get_var(fid, varid, CF_LAT)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, ISTAT)

        WRITE(WINDBG%FHNDL,*) 'MIN,MAX CF_LON=', MINVAL(CF_LON), MAXVAL(CF_LON)
        WRITE(WINDBG%FHNDL,*) 'MIN,MAX CF_LAT=', MINVAL(CF_LAT), MAXVAL(CF_LAT)
        FLUSH(WINDBG%FHNDL)

        TheInfo % nx_dim = NDX_WIND_FD
        TheInfo % ny_dim = NDY_WIND_FD
        allocate(TheInfo % LON(NDX_WIND_FD, NDY_WIND_FD), TheInfo % LAT(NDX_WIND_FD, NDY_WIND_FD), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 47')
        DO IX=1,NDX_WIND_FD
          DO IY=1,NDY_WIND_FD
            TheInfo % LON(IX,IY) = CF_LON(IX,IY)
            TheInfo % LAT(IX,IY) = CF_LAT(IX,IY)
          END DO
        END DO
        DEALLOCATE(CF_LON, CF_LAT)
        CALL COMPUTE_CF_COEFFICIENTS(TheInfo)
        Deallocate(TheInfo % LON, TheInfo % LAT)
# ifdef MPI_PARALL_GRID
      END IF
# endif
      END SUBROUTINE
!**********************************************************************
!*    This is for direct to elements forcing in netcdf                *
!*    wind_time  as usual                                             *
!*    float Uwind(wind_time, mnp)                                     *
!*    float Vwind(wind_time, mnp)                                     *
!*    should be better than text.                                     *
!**********************************************************************
      SUBROUTINE READ_DIRECT_NETCDF_CF(eVAR, RECORD_IN, outwind)
      USE NETCDF
      USE DATAPOOL
      IMPLICIT NONE
      TYPE(VAR_NETCDF_CF)                :: eVAR
      INTEGER, INTENT(in)                :: RECORD_IN
      REAL(rkind), INTENT(out)           :: outwind(MNP,2)
      character (len = *), parameter :: CallFct="READ_DIRECT_NETCDF_CF"
      INTEGER                            :: FID, ID
      real(rkind) :: U_tot(np_total), V_tot(np_total)
      real(rkind) :: Vtotal1(np_total), Vtotal2(np_total)
      real(rkind) :: Vlocal(MNP)
      real(rkind) :: cf_scale_factor, cf_add_offset
# ifdef MPI_PARALL_GRID
      integer IP_glob, IP
# endif
      character(len=10) :: eStrU, eStrV
      cf_scale_factor = eVAR % cf_scale_factor
      cf_add_offset = eVAR % cf_add_offset
      IF (eVAR % idVar == 1) THEN
        eStrU='Uwind'
        eStrV='Vwind'
      ELSE
        eStrU='UsurfCurr'
        eStrV='VsurfCurr'
      END IF
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
# endif
        CALL TEST_FILE_EXIST_DIE("Missing file : ", TRIM(eVAR % eFileName))
        ISTAT = NF90_OPEN(TRIM(eVAR % eFileName), NF90_NOWRITE, FID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)

        ISTAT = NF90_inq_varid(FID, TRIM(eStrU), ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)

        ISTAT = NF90_GET_VAR(FID, ID, U_tot, start = (/ 1, RECORD_IN /), count = (/ np_total, 1 /))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, ISTAT)

        ISTAT = NF90_inq_varid(FID, TRIM(eStrV), ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, ISTAT)

        ISTAT = NF90_GET_VAR(FID, ID, V_tot, start = (/ 1, RECORD_IN /), count = (/ np_total, 1 /))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)

        ISTAT = NF90_CLOSE(FID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, ISTAT)

        Vtotal1 = cf_add_offset + cf_scale_factor*U_tot
        Vtotal2 = cf_add_offset + cf_scale_factor*V_tot
# ifdef MPI_PARALL_GRID
      END IF
# endif
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND) THEN
        DO IP=1,MNP
          IP_glob=iplg(IP)
          outwind(IP,1)=Vtotal1(IP_glob)
          outwind(IP,2)=Vtotal2(IP_glob)
        END DO
      ELSE
        CALL SCATTER_ONED_ARRAY(Vtotal1, Vlocal)
        outwind(:,1)=Vlocal
        CALL SCATTER_ONED_ARRAY(Vtotal2, Vlocal)
        outwind(:,2)=Vlocal
      END IF
# else
      outwind(:,1)=Vtotal1
      outwind(:,2)=Vtotal2
# endif
      WRITE(WINDBG%FHNDL,*) 'READ_DIRECT_NETCDF_CF'
      WRITE(WINDBG%FHNDL,*) 'RECORD_IN=', RECORD_IN
      WRITE(WINDBG%FHNDL,*) 'UWIND_FE, min/max=', minval(outwind(:,1)), maxval(outwind(:,1))
      WRITE(WINDBG%FHNDL,*) 'VWIND_FE, min/max=', minval(outwind(:,2)), maxval(outwind(:,2))
      FLUSH(WINDBG%FHNDL)
      END SUBROUTINE
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE READ_DIRECT_NETCDF_CF1(eVAR, RECORD_IN, outvar)
      USE NETCDF
      USE DATAPOOL
      IMPLICIT NONE
      TYPE(VAR_NETCDF_CF)                :: eVAR
      INTEGER, INTENT(in)                :: RECORD_IN
      REAL(rkind), INTENT(out)           :: outvar(MNP)
      character (len = *), parameter :: CallFct="READ_DIRECT_NETCDF_CF"
      INTEGER                            :: FID, ID
      real(rkind) :: VAR_tot(np_total)
      real(rkind) :: Vtotal(np_total)
      real(rkind) :: Vlocal(MNP)
      real(rkind) :: cf_scale_factor, cf_add_offset
# ifdef MPI_PARALL_GRID
      integer IP_glob, IP
# endif
      character(len=10) :: eStrU, eStrV
      cf_scale_factor = eVAR % cf_scale_factor
      cf_add_offset = eVAR % cf_add_offset
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
# endif
        CALL TEST_FILE_EXIST_DIE("Missing file : ", TRIM(eVAR % eFileName))
        ISTAT = NF90_OPEN(TRIM(eVAR % eFileName), NF90_NOWRITE, FID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)

        ISTAT = NF90_inq_varid(FID, TRIM(eVAR % eString), ID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)

        ISTAT = NF90_GET_VAR(FID, ID, VAR_tot, start = (/ 1, RECORD_IN /), count = (/ np_total, 1 /))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, ISTAT)

        ISTAT = NF90_CLOSE(FID)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, ISTAT)

        Vtotal = cf_add_offset + cf_scale_factor*VAR_tot
# ifdef MPI_PARALL_GRID
      END IF
# endif
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND) THEN
        DO IP=1,MNP
          IP_glob=iplg(IP)
          outvar(IP)=Vtotal(IP_glob)
        END DO
      ELSE
        CALL SCATTER_ONED_ARRAY(Vtotal, outvar)
      END IF
# else
      outvar=Vtotal
# endif
      WRITE(WINDBG%FHNDL,*) 'READ_DIRECT_NETCDF_CF1'
      WRITE(WINDBG%FHNDL,*) 'RECORD_IN=', RECORD_IN
      FLUSH(WINDBG%FHNDL)
      END SUBROUTINE
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE INIT_DIRECT_NETCDF_CF(eVAR, MULTIPLE_IN, eFileName, eString)
      USE NETCDF
      USE DATAPOOL
      IMPLICIT NONE
      TYPE(VAR_NETCDF_CF), intent(inout) :: eVAR
      logical, intent(in) :: MULTIPLE_IN
      character(len=100), intent(in) :: eFileName
      character(len=*), intent(in) :: eString
!      
      INTEGER           :: fid, varid
      INTEGER           :: dimidsB(3), dimids(3)
      integer nbChar
      character (len=20) :: WindTimeStr
      character(len=100) :: CHRERR
      character (len = *), parameter :: CallFct="INIT_DIRECT_NETCDF_CF"
      character (len=100) :: eStrUnitTime
      real(rkind) :: ConvertToDay
      real(rkind) :: eTimeStart
      real(rkind) :: cf_scale_factor, cf_add_offset
      integer nbtime_mjd
      real(rkind), allocatable :: wind_time_mjd(:)
      integer eInt(1)
      real(rkind) :: eReal(2)
      integer IPROC
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN .or. (myrank .eq. 0)) THEN
# endif
        ISTAT = nf90_open(TRIM(eFileName), nf90_nowrite, fid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)

        ! Reading wind attributes

!        Print *, 'eString=', TRIM(eString)
!        Print *, 'FNAME=', TRIM(eFileName)
        ISTAT = nf90_inq_varid(fid, TRIM(eString), varid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)

        ISTAT = nf90_inquire_variable(fid, varid, dimids=dimidsB)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, ISTAT)

        ISTAT = nf90_inquire_dimension(fid, dimidsB(3), name=WindTimeStr)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, ISTAT)
        WRITE(WINDBG%FHNDL,*) 'variable used for time=', TRIM(WindTimeStr)
        FLUSH(WINDBG%FHNDL)

        ISTAT = nf90_get_att(fid, varid, "scale_factor", cf_scale_factor)
        IF (ISTAT /= 0) THEN
          CHRERR = nf90_strerror(ISTAT)
          WRITE(WINDBG%FHNDL,*) 'CHRERR=', TRIM(CHRERR)
          cf_scale_factor=ONE
        ENDIF
        WRITE(WINDBG%FHNDL,*) 'cf_scale_factor=', cf_scale_factor
        FLUSH(WINDBG%FHNDL)

        ISTAT = nf90_get_att(fid, varid, "add_offset", cf_add_offset)
        IF (ISTAT /= 0) THEN
          CHRERR = nf90_strerror(ISTAT)
          WRITE(WINDBG%FHNDL,*) 'CHRERR=', TRIM(CHRERR)
          cf_add_offset=ZERO
        ENDIF
        WRITE(WINDBG%FHNDL,*) 'cf_add_offset=', cf_add_offset
        FLUSH(WINDBG%FHNDL)

        ! Reading time
       
        ISTAT = nf90_inq_varid(fid, WindTimeStr, varid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)

        ISTAT = nf90_inquire_attribute(fid, varid, "units", len=nbChar)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, ISTAT)

        ISTAT = nf90_get_att(fid, varid, "units", eStrUnitTime)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, ISTAT)
        CALL CF_EXTRACT_TIME(eStrUnitTime, ConvertToDay, eTimeStart)
        WRITE(WINDBG%FHNDL,*) 'eTimeStart=', eTimeStart
        FLUSH(WINDBG%FHNDL)

        ISTAT = nf90_inquire_variable(fid, varid, dimids=dimids)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, ISTAT)

        ISTAT = nf90_inquire_dimension(fid, dimids(1), len=nbtime_mjd)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, ISTAT)

        allocate(wind_time_mjd(nbtime_mjd), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 48')

        ISTAT = nf90_get_var(fid, varid, wind_time_mjd)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, ISTAT)

        ISTAT = nf90_close(fid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)

        wind_time_mjd(:) = wind_time_mjd(:)*ConvertToDay + eTimeStart
        CALL CHECK_WIND_TIME(nbtime_mjd, WIND_TIME_MJD)
# ifdef MPI_PARALL_GRID
      END IF
# endif
#ifdef MPI_PARALL_GRID
      IF (.NOT. MULTIPLE_IN_WIND) THEN
        IF (myrank .eq. 0) THEN
          eInt(1)=nbtime_mjd
          DO IPROC=2,nproc
            CALL MPI_SEND(eInt,1,itype, iProc-1, 811, comm, ierr)
          END DO
          eReal(1)=cf_scale_factor
          eReal(2)=cf_add_offset
          DO IPROC=2,nproc
            CALL MPI_SEND(eReal,2,rtype, iProc-1, 813, comm, ierr)
          END DO
          DO IPROC=2,nproc
            CALL MPI_SEND(WIND_TIME_MJD,nbtime_mjd,rtype, iProc-1, 812, comm, ierr)
          END DO
        ELSE
          CALL MPI_RECV(eInt,1,itype, 0, 811, comm, istatus, ierr)
          nbtime_mjd=eInt(1)
          CALL MPI_RECV(eReal,2,rtype, 0, 813, comm, istatus, ierr)
          cf_scale_factor=eReal(1)
          cf_add_offset=eReal(2)
          ALLOCATE(WIND_TIME_MJD(nbtime_mjd), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 3')
          CALL MPI_RECV(WIND_TIME_MJD,nbtime_mjd,rtype, 0, 812, comm, istatus, ierr)
        END IF
      END IF
#endif
      eVAR % cf_scale_factor = cf_scale_factor
      eVAR % cf_add_offset = cf_add_offset
      eVAR % nbTime = nbtime_mjd
      allocate(eVAR % ListTime(nbtime_mjd))
      eVAR % ListTime = WIND_TIME_MJD
      eVAR % eFileName = eFileName
      eVAR % eString = eString
      IF (TRIM(eString) == 'Uwind') THEN
        eVAR % idVar = 1
      ELSE
        eVAR % idVar = 2
      END IF

      END SUBROUTINE
#endif
#ifdef GRIB_API_ECMWF
!****************************************************************************
!* Raw reading of time entry for GRIB                                       *
!****************************************************************************
      SUBROUTINE RAW_READ_TIME_OF_GRIB_FILE(ifile, eGrib, STEPRANGE_IN, eTimeOut)
      USE DATAPOOL
      USE GRIB_API
      IMPLICIT NONE
      integer, intent(in) :: ifile
      integer, intent(in) :: eGrib
      LOGICAL, intent(in) :: STEPRANGE_IN
      real(rkind), intent(out) :: eTimeOut
      !
      LOGICAL :: USE_DATATIME = .TRUE.
      integer eYear, eMonth, eDay, resYear, resMonth
      integer eHour, eMin, eSec, resHour, resMin
      integer dataDate, stepRange, dataTime
      character (len=15) :: eStrTime
      REAL(rkind) :: eTimeBase
      call grib_get(eGrib, 'dataDate', dataDate)
      eYear=(dataDate - mod(dataDate,10000))/10000
      resYear=dataDate - 10000*eYear
      eMonth=(resYear - mod(resYear,100))/100
      resMonth=resYear - 100*eMonth;
      eDay=resMonth
      IF (STEPRANGE_IN) THEN
        call grib_get(eGrib, 'stepRange', stepRange)
      ELSE
        stepRange=0
      END IF
      IF (USE_DATATIME) THEN
        call grib_get(eGrib, 'dataTime', dataTime)
        eHour=(dataTime - mod(dataTime,100))/100
        eMin=dataTime - 100*eHour
        eSec=0
      ELSE
        eHour=0
        eMin=0
        eSec=0
      END IF
      WRITE(eStrTime,10) eYear, eMonth, eDay, eHour, eMin, eSec
 10   FORMAT(i4.4,i2.2,i2.2,'.',i2.2,i2.2,i2.2)
      CALL CT2MJD(eStrTime, eTimeBase)
      eTimeOut=eTimeBase + MyREAL(stepRange)/24.0_rkind
      END SUBROUTINE
!****************************************************************************
!* Reading time from a GRIB file                                            *
!****************************************************************************
      SUBROUTINE READ_TIME_OF_GRIB_FILE(eTimeOut, eFile, STEPRANGE_IN)
      USE DATAPOOL
      USE GRIB_API
      IMPLICIT NONE
      REAL(rkind), intent(out) :: eTimeOut
      CHARACTER(len=140), intent(in) :: eFile
      LOGICAL, intent(in) :: STEPRANGE_IN
      integer ifile, i, n
      integer, allocatable :: igrib(:)
      CALL TEST_FILE_EXIST_DIE("Missing grib file: ", TRIM(eFile))
      CALL GRIB_OPEN_FILE(ifile, TRIM(eFile), 'r')
      call grib_count_in_file(ifile,n)
      allocate(igrib(n))
      i=1
      call grib_new_from_file(ifile, igrib(i))
      CALL RAW_READ_TIME_OF_GRIB_FILE(ifile, igrib(i), STEPRANGE_IN, eTimeOut)
      CALL GRIB_CLOSE_FILE(ifile)
      deallocate(igrib)
      END SUBROUTINE
!****************************************************************************
!* Reading grid information from a GRIB file (case 1)                       *
!****************************************************************************
      SUBROUTINE READ_GRID_INFO_FROM_GRIB_TYPE1(TheInfo, eGrib)
      USE DATAPOOL
      USE GRIB_API
      IMPLICIT NONE
      type(FD_FORCING_GRID), intent(out) :: TheInfo
      integer, intent(in) :: eGrib
      !
      REAL(rkind) ::longitudeOfFirstPointInDegrees, latitudeOfFirstPointInDegrees, longitudeOfLastPointInDegrees, latitudeOfLastPointInDegrees
      REAL(rkind) :: deltaLAT, deltaLON
      REAL(rkind) :: iDirectionIncrement, jDirectionIncrement
      integer nx_dim, ny_dim, iX, iY
      call grib_get(eGrib,"numberOfPointsAlongAParallel", nx_dim)
      call grib_get(eGrib,"numberOfPointsAlongAMeridian", ny_dim)
      WRITE(STAT%FHNDL, *) 'nx_dim=', nx_dim
      WRITE(STAT%FHNDL, *) 'ny_dim=', ny_dim
      TheInfo % nx_dim = nx_dim
      TheInfo % ny_dim = ny_dim
      allocate(TheInfo % LON(nx_dim, ny_dim), TheInfo % LAT(nx_dim, ny_dim), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 47')
      call grib_get(eGrib, 'longitudeOfFirstGridPointInDegrees', longitudeOfFirstPointInDegrees)
      call grib_get(eGrib, 'latitudeOfFirstGridPointInDegrees', latitudeOfFirstPointInDegrees)
      call grib_get(eGrib, 'longitudeOfLastGridPointInDegrees', longitudeOfLastPointInDegrees)
      call grib_get(eGrib, 'latitudeOfLastGridPointInDegrees', latitudeOfLastPointInDegrees)

      call grib_get(eGrib, 'iDirectionIncrementInDegrees', iDirectionIncrement)
      call grib_get(eGrib, 'jDirectionIncrementInDegrees', jDirectionIncrement)

      WRITE(STAT%FHNDL, *) 'LONGITUDE'
      WRITE(STAT%FHNDL, *) 'longitudeOfFirstGridPointInDegrees=', longitudeOfFirstPointInDegrees
      WRITE(STAT%FHNDL, *) 'longitudeOfLastGridPointInDegrees=', longitudeOfLastPointInDegrees
      WRITE(STAT%FHNDL, *) 'LATITUDE'
      WRITE(STAT%FHNDL, *) 'latitudeOfFirstGridPointInDegrees=', latitudeOfFirstPointInDegrees
      WRITE(STAT%FHNDL, *) 'latitudeOfLastGridPointInDegrees=', latitudeOfLastPointInDegrees

      WRITE(STAT%FHNDL, *) 'iDirectionIncrement=', iDirectionIncrement
      WRITE(STAT%FHNDL, *) 'jDirectionIncrement=', jDirectionIncrement
      deltaLON=(longitudeOfLastPointInDegrees - longitudeOfFirstPointInDegrees)/(nx_dim - 1)
      deltaLAT=(latitudeOfLastPointInDegrees - latitudeOfFirstPointInDegrees)/(ny_dim - 1)
      DO iX=1,nx_dim
        DO iY=1,ny_dim
          TheInfo % LON(iX,iY)=longitudeOfFirstPointInDegrees + (iX-1)*deltaLON
          TheInfo % LAT(iX,iY)=latitudeOfFirstPointInDegrees + (iY-1)*deltaLAT
        END DO
      END DO
      END SUBROUTINE
!****************************************************************************
!* Reading grid information from a GRIB file (case 2)                       *
!****************************************************************************
      SUBROUTINE phirot2phi(eLatOut, phirot, rlarot, polphi, pollam, polgam)
      USE DATAPOOL
      IMPLICIT NONE
      real(rkind), intent(out) :: eLatOut
      real(rkind), intent(in) :: phirot, rlarot, polphi, pollam, polgam
      real(rkind) :: zsinpol
      real(rkind) :: zcospol
      real(rkind) :: zphis
      real(rkind) :: zrlas
      real(rkind) :: zarg
      real(rkind) :: zgam
      zsinpol = sin(DEGRAD * polphi)
      zcospol = cos(DEGRAD * polphi)
      zphis  = DEGRAD * phirot
      
      if (rlarot .gt. 180.) THEN
        zrlas = rlarot - 360.
      ELSE
        zrlas = rlarot
      END IF
      zrlas = DEGRAD * zrlas;
      if (ABS(polgam) .gt. 0) THEN
        zgam  = DEGRAD * polgam;
        zarg = zsinpol*sin(zphis) + zcospol*cos(zphis) * ( cos(zrlas)*cos(zgam) - sin(zgam)*sin(zrlas))
      ELSE
        zarg = zcospol * cos(zphis) * cos(zrlas) + zsinpol * sin(zphis)
      END IF
      eLatOut = RADDEG * asin(zarg);
      END SUBROUTINE
!****************************************************************************
!* Reading grid information from a GRIB file (case 2)                       *
!****************************************************************************
      SUBROUTINE rlarot2rla(eLonOut, phirot, rlarot, polphi, pollam, polgam)
      USE DATAPOOL
      IMPLICIT NONE
      real(rkind), intent(out) :: eLonOut
      real(rkind), intent(in) :: phirot, rlarot, polphi, pollam, polgam
      !
      real(rkind) :: zsinpol
      real(rkind) :: zcospol
      real(rkind) :: zphis
      real(rkind) :: zrlas
      real(rkind) :: zgam, zlampol
      real(rkind) :: zarg1, zarg2
      zsinpol = sin(DEGRAD * polphi)
      zcospol = cos(DEGRAD * polphi)
      zphis  = DEGRAD * phirot

      
      if (rlarot .gt. 180.) THEN
        zrlas = rlarot - 360.
      ELSE
        zrlas = rlarot
      END IF
      zrlas = DEGRAD * zrlas
      zlampol = DEGRAD * pollam
      IF (ABS(polgam) .gt. 0) THEN
        zgam    = DEGRAD * polgam;
        zarg1   = sin (zlampol) *                                                       &
   &     (- zsinpol*cos(zphis) * (cos(zrlas)*cos(zgam) - sin(zrlas)*sin(zgam))          &
   &     + zcospol * sin(zphis))                                                        &
   &     - cos (zlampol)*cos(zphis) * (sin(zrlas)*cos(zgam) + cos(zrlas)*sin(zgam))
        zarg2   = cos (zlampol) *                                                       &
   &     (- zsinpol*cos(zphis) * (cos(zrlas)*cos(zgam) - sin(zrlas)*sin(zgam))          &
   &     + zcospol * sin(zphis))                                                        &
   &     + sin (zlampol)*cos(zphis) * (sin(zrlas)*cos(zgam) + cos(zrlas)*sin(zgam))
      ELSE
        zarg1   = sin (zlampol) * (-zsinpol * cos(zrlas) * cos(zphis)  +                &
   &     zcospol *              sin(zphis)) -                                           &
   &     cos (zlampol) *             sin(zrlas) * cos(zphis)
        zarg2   = cos (zlampol) * (-zsinpol * cos(zrlas) * cos(zphis)  +                &
   &     zcospol *              sin(zphis)) +                                           &
   &     sin (zlampol) *             sin(zrlas) * cos(zphis)
      END IF
      if (zarg2 .eq. 0) zarg2=1.0e-20
      eLonOut = RADDEG * atan2(zarg1,zarg2);
      END SUBROUTINE
!****************************************************************************
!* Reading grid information from a GRIB file (case 2)                       *
!****************************************************************************
      SUBROUTINE READ_GRID_INFO_FROM_GRIB_TYPE2(TheInfo, eGrib)
      USE DATAPOOL
      USE GRIB_API
      IMPLICIT NONE
      type(FD_FORCING_GRID), intent(out) :: TheInfo
      integer, intent(in) :: eGrib
      real(rkind) :: latitudeOfSouthernPoleInDegrees, longitudeOfSouthernPoleInDegrees
      real(rkind) :: angleOfRotationInDegrees
      real(rkind) :: latitudeOfFirstGridPointInDegrees, longitudeOfFirstGridPointInDegrees
      real(rkind) :: latitudeOfLastGridPointInDegrees, longitudeOfLastGridPointInDegrees
      real(rkind) :: iDirectionIncrementInDegrees, jDirectionIncrementInDegrees
      real(rkind) :: pollat_sp, pollon_sp, polgam, zstartlon_tot, zstartlat_tot
      real(rkind) :: zendlon_tot, zendlat_tot, dlon, dlat
      real(rkind) :: eLonR, eLatR, eLonOut, eLatOut
      real(rkind) :: pollat, pollon
      real(rkind) :: startlon_tot, startlat_tot
      integer :: nx_dim, ny_dim, iX, iY
      !
      call grib_get(eGrib, 'latitudeOfSouthernPoleInDegrees',latitudeOfSouthernPoleInDegrees)
      call grib_get(eGrib, 'longitudeOfSouthernPoleInDegrees',longitudeOfSouthernPoleInDegrees)
      call grib_get(eGrib, 'angleOfRotationInDegrees',angleOfRotationInDegrees)
      call grib_get(eGrib, 'latitudeOfFirstGridPointInDegrees',latitudeOfFirstGridPointInDegrees)
      call grib_get(eGrib, 'longitudeOfFirstGridPointInDegrees',longitudeOfFirstGridPointInDegrees)
      call grib_get(eGrib, 'latitudeOfLastGridPointInDegrees',latitudeOfLastGridPointInDegrees)
      call grib_get(eGrib, 'longitudeOfLastGridPointInDegrees',longitudeOfLastGridPointInDegrees)
      call grib_get(eGrib, 'iDirectionIncrementInDegrees',iDirectionIncrementInDegrees)
      call grib_get(eGrib, 'jDirectionIncrementInDegrees',jDirectionIncrementInDegrees)
      pollat_sp=latitudeOfSouthernPoleInDegrees
      pollon_sp=longitudeOfSouthernPoleInDegrees
      polgam=angleOfRotationInDegrees
      zstartlon_tot=longitudeOfFirstGridPointInDegrees
      zstartlat_tot=latitudeOfFirstGridPointInDegrees
      zendlon_tot=longitudeOfLastGridPointInDegrees
      zendlat_tot=latitudeOfLastGridPointInDegrees
      dlon=iDirectionIncrementInDegrees
      dlat=jDirectionIncrementInDegrees
      !
      ! Now reading the mapped grid
      !
      call grib_get(eGrib,"Nx", nx_dim)
      call grib_get(eGrib,"Ny", ny_dim)
      WRITE(STAT%FHNDL, *) 'nx_dim = ', nx_dim
      WRITE(STAT%FHNDL, *) 'ny_dim = ', ny_dim
      TheInfo % nx_dim = nx_dim
      TheInfo % ny_dim = ny_dim
      pollat= - pollat_sp
      pollon= pollon_sp - 180.
      startlon_tot=zstartlon_tot
      startlat_tot=zstartlat_tot
      
      allocate(TheInfo % LON(nx_dim, ny_dim), TheInfo % LAT(nx_dim, ny_dim), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 47')
      DO iX=1,nx_dim
        DO iY=1,ny_dim
          eLonR = startlon_tot + MyREAL(iX-1)*dlon
          eLatR = startlat_tot + MyREAL(iY-1)*dlat
          CALL phirot2phi(eLatOut, eLatR, eLonR, pollat, pollon, polgam)
          CALL rlarot2rla(eLonOut, eLatR, eLonR, pollat, pollon, polgam)
          TheInfo % LON(iX,iY) = eLonOut
          TheInfo % LAT(iX,iY) = eLatOut
        END DO
      END DO
      END SUBROUTINE
!****************************************************************************
!* Reading grid information from a GRIB file (case 1)                       *
!****************************************************************************
      SUBROUTINE READ_GRID_INFO_FROM_GRIB_TYPE3(TheInfo, eGrib)
      USE DATAPOOL
      USE GRIB_API
      IMPLICIT NONE
      type(FD_FORCING_GRID), intent(out) :: TheInfo
      integer, intent(in) :: eGrib
      REAL(rkind), allocatable :: LON_serial(:), LAT_serial(:), DATA_Serial(:)
      integer nx_dim, ny_dim, idx, eProd
      integer status, iX, iY
      !
      call grib_get(eGrib,"Nx", nx_dim)
      call grib_get(eGrib,"Ny", ny_dim)
      WRITE(STAT%FHNDL, *) 'nx_dim = ', nx_dim
      WRITE(STAT%FHNDL, *) 'ny_dim = ', ny_dim
      TheInfo % nx_dim = nx_dim
      TheInfo % ny_dim = ny_dim
      allocate(TheInfo % LON(nx_dim, ny_dim), TheInfo % LAT(nx_dim, ny_dim), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 47')
      eProd=nx_dim*ny_dim
      allocate(LON_serial(eProd), LAT_serial(eProd), DATA_serial(eProd), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 48')
      call grib_get_data(eGrib, LAT_serial, LON_serial, DATA_serial, status)
      idx=0
      DO iY=1,ny_dim
        DO iX=1,nx_dim
          idx=idx+1
          TheInfo % LON(iX,iY)=LON_serial(idx)
          TheInfo % LAT(iX,iY)=LAT_serial(idx)
        END DO
      END DO
      DEALLOCATE(LON_serial, LAT_serial, DATA_serial)      
      END SUBROUTINE
!****************************************************************************
!* Reading grid information from a GRIB file                                *
!****************************************************************************
      SUBROUTINE READ_GRID_INFO_FROM_GRIB(TheInfo, TheFile, shortName, GRIB_TYPE)
      USE DATAPOOL
      USE GRIB_API
      IMPLICIT NONE
      type(FD_FORCING_GRID), intent(out) :: TheInfo
      character(len=140), intent(in) :: TheFile
      character(len=20), intent(in) :: shortName
      integer, intent(in) :: GRIB_TYPE
      !
      integer ifile, i, n
      logical WeFound
      integer, allocatable :: igrib(:)
      character(len=100) eShortName
      integer eProd
      integer status, idx
      integer iX, iY
      integer nx_dim, ny_dim
      !
      CALL TEST_FILE_EXIST_DIE("Missing grib file: ", TRIM(TheFile))
      CALL GRIB_OPEN_FILE(ifile, TRIM(TheFile), 'r')
      call grib_count_in_file(ifile,n)
      allocate(igrib(n))
      !
      WeFound=.FALSE.;
      DO i=1,n
        call grib_new_from_file(ifile, igrib(i))
        call grib_get(igrib(i), 'shortName', eShortName)
        IF ((TRIM(eShortName) .eq. shortName).and.(WeFound .eqv. .FALSE.)) THEN
          IF (GRIB_FILE_TYPE .eq. 1) THEN
            CALL READ_GRID_INFO_FROM_GRIB_TYPE1(TheInfo, igrib(i))
          END IF
          IF (GRIB_FILE_TYPE .eq. 2) THEN
            CALL READ_GRID_INFO_FROM_GRIB_TYPE2(TheInfo, igrib(i))
          END IF
          IF (GRIB_FILE_TYPE .eq. 3) THEN
            CALL READ_GRID_INFO_FROM_GRIB_TYPE3(TheInfo, igrib(i))
          END IF
          WeFound=.TRUE.
        END IF
      END DO
      IF (WeFound .eqv. .FALSE.) THEN
        Print *, 'Failed to find the wind variable in the grib file'
        CALL WWM_ABORT("Wind has not been found in grib file")          
      END IF
      WRITE(STAT%FHNDL, *) 'WeFound=', WeFound
      CALL GRIB_CLOSE_FILE(ifile)
      deallocate(igrib)
      END SUBROUTINE
!****************************************************************************
!* This is functionality for reading GRIB file wind input                   *
!* Specific supported cases (or wished): ECMWF (IFS), COSMO, DHMZ (ALADIN)  *
!****************************************************************************
      SUBROUTINE INIT_GRIB_WIND
      USE DATAPOOL
      USE GRIB_API
      IMPLICIT NONE
      INTEGER IT
      INTEGER ifile, i, n
      integer, allocatable :: igrib(:)
      integer WeFound
      REAL(rkind) :: eTimeMjd
      integer IPROC, eInt(1)
      integer iX, iY
      integer nbtime_mjd
      integer eProd
      character(len=20) shortName
      REAL(rkind), allocatable :: wind_time_mjd(:)
      REAL cf_scale_factor, cf_add_offset
      TYPE(FD_FORCING_GRID) :: TheInfo
      integer GRIB_TYPE
!     PRint *, 'Begin of INIT_GRIB_WIND'
      WRITE(WINDBG%FHNDL, *) 'GRIB_FILE_TYPE=', GRIB_FILE_TYPE
      WRITE(WINDBG%FHNDL, *) 'MULTIPLE_IN_WIND=', MULTIPLE_IN_WIND
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
# endif
        OPEN(WIN%FHNDL,FILE=WIN%FNAME,STATUS='OLD',ACTION='READ', IOSTAT = ISTAT)
        NUM_GRIB_FILES = 0
        DO
          READ( WIN%FHNDL, *, IOSTAT = ISTAT )
          IF ( ISTAT /= 0 ) EXIT
          NUM_GRIB_FILES = NUM_GRIB_FILES + 1
        END DO
        WRITE(WINDBG%FHNDL,*) 'NUM_GRIB_FILES=', NUM_GRIB_FILES
        REWIND (WIN%FHNDL)
        !
        ALLOCATE(GRIB_FILE_NAMES(NUM_GRIB_FILES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 18')
        DO IT = 1, NUM_GRIB_FILES
          READ(WIN%FHNDL, *) GRIB_FILE_NAMES(IT)
          WRITE(WINDBG%FHNDL,*) IT, TRIM(GRIB_FILE_NAMES(IT))
        END DO
        CLOSE (WIN%FHNDL)
        FLUSH(WINDBG%FHNDL)
        !
        nbtime_mjd=NUM_GRIB_FILES
        allocate(wind_time_mjd(nbtime_mjd), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 48')
        DO IT=1, nbTime_mjd
          WRITE(WINDBG%FHNDL, *) '---------------------------------------'
          WRITE(WINDBG%FHNDL, *) 'IT=', IT, 'file = ',  TRIM(GRIB_FILE_NAMES(IT))
          CALL READ_TIME_OF_GRIB_FILE(eTimeMjd, GRIB_FILE_NAMES(IT), USE_STEPRANGE)
          wind_time_mjd(IT)=eTimeMjd
        END DO
        FLUSH(WINDBG%FHNDL)
        cf_scale_factor=ONE
        cf_add_offset=ZERO
        !
        ! Now the longitude/latitude to read.
        !
        shortName='10u'
        GRIB_TYPE = GRIB_FILE_TYPE
        IT=1
        CALL READ_GRID_INFO_FROM_GRIB(TheInfo, GRIB_FILE_NAMES(IT), shortName, GRIB_TYPE)
        NDX_WIND_FD = TheInfo % nx_dim
        NDY_WIND_FD = TheInfo % ny_dim
        allocate(UWIND_FD(NDX_WIND_FD, NDY_WIND_FD), VWIND_FD(NDX_WIND_FD, NDY_WIND_FD), stat=istat)
        CALL COMPUTE_CF_COEFFICIENTS(TheInfo)
        DEALLOCATE(TheInfo % LON, TheInfo % LAT)
# ifdef MPI_PARALL_GRID
      END IF
# endif
# ifdef MPI_PARALL_GRID
      IF (.NOT. MULTIPLE_IN_WIND) THEN
        IF (myrank .eq. 0) THEN
          eInt(1)=nbtime_mjd
          DO IPROC=2,nproc
            CALL MPI_SEND(eInt,1,itype, iProc-1, 811, comm, ierr)
          END DO
          DO IPROC=2,nproc
            CALL MPI_SEND(WIND_TIME_MJD,nbtime_mjd,rtype, iProc-1, 812, comm, ierr)
          END DO
        ELSE
          CALL MPI_RECV(eInt,1,itype, 0, 811, comm, istatus, ierr)
          nbtime_mjd=eInt(1)
          ALLOCATE(WIND_TIME_MJD(nbtime_mjd), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 3')
          CALL MPI_RECV(WIND_TIME_MJD,nbtime_mjd,rtype, 0, 812, comm, istatus, ierr)
        END IF
      END IF
# endif
      eVAR_WIND % cf_scale_factor = cf_scale_factor
      eVAR_WIND % cf_add_offset = cf_add_offset
      eVAR_WIND % nbTime= nbtime_mjd
      allocate(eVAR_WIND % ListTime(nbtime_mjd))
      eVAR_WIND % ListTime = WIND_TIME_MJD
!     PRint *, 'End of INIT_GRIB_WIND'
      END SUBROUTINE
!****************************************************************************
!* The read subroutine                                                      *
!****************************************************************************
      SUBROUTINE READ_GRIB_WIND(IT, outwind)
      USE DATAPOOL
      USE GRIB_API
      IMPLICIT NONE
      integer, intent(in) :: IT
      REAL(rkind), INTENT(out)           :: outwind(MNP,2)
      REAL(rkind)                        :: outTotal(MNP_WIND,2)
      REAL(rkind)                        :: Vtotal(np_total)
      REAL(rkind)                        :: Vlocal(MNP)
      INTEGER ifile, irec, n, iret
      integer, allocatable :: igrib(:)
      integer WeFoundU, WeFoundV
      integer i, j, idx
      character(len=100) eShortName
      real(rkind) valueU(NDX_WIND_FD*NDY_WIND_FD)
      real(rkind) valueV(NDX_WIND_FD*NDY_WIND_FD)
!     Print *, 'NDX_WIND_FD=', NDX_WIND_FD
!     Print *, 'NDY_WIND_FD=', NDY_WIND_FD
      !
!     Print *, 'Begin of READ_GRIB_WIND'
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
# endif
        WRITE(WINDBG%FHNDL,*) 'IT=', IT, 'file = ',  TRIM(GRIB_FILE_NAMES(IT))
!       Print *, 'GRIB_FILE_NAMES(IT)=', TRIM(GRIB_FILE_NAMES(IT))
        CALL TEST_FILE_EXIST_DIE("Missing grib file: ", TRIM(GRIB_FILE_NAMES(IT)))
        CALL GRIB_OPEN_FILE(ifile, TRIM(GRIB_FILE_NAMES(IT)), 'r')
        call grib_count_in_file(ifile,n)
        WRITE(WINDBG%FHNDL,*) 'n=', n
        allocate(igrib(n))
        WeFoundU=0
        WeFoundV=0
        DO irec=1,n
          call grib_new_from_file(ifile, igrib(irec), iret)
!         Print *, 'Step 1'
          call grib_get(igrib(irec), 'shortName', eShortName)
!         Print *, 'Step 2'
!         Print *, 'eShortName=', eShortName
          IF ((TRIM(eShortName) .eq. '10u').and.(WeFoundU .eq. 0)) THEN
            WeFoundU=1
!           Print *, 'Step 3'
            CALL grib_get(igrib(irec), 'values', valueU)
!           Print *, 'Step 4'
          END IF
          IF ((TRIM(eShortName) .eq. '10v').and.(WeFoundV .eq. 0)) THEN
            WeFoundV=1
!           Print *, 'Step 5'
            CALL grib_get(igrib(irec), 'values', valueV)
!           Print *, 'Step 6'
          END IF
        END DO
        idx=0
        DO J=1,NDY_WIND_FD
          DO I=1,NDX_WIND_FD
            idx=idx+1
            UWIND_FD(I,J)=valueU(idx)
            VWIND_FD(I,J)=valueV(idx)
          END DO
        END DO
        CALL GRIB_CLOSE_FILE(ifile)
        CALL KERNEL_INTERP_UV_WINDFD(outTotal)
# ifdef MPI_PARALL_GRID
      END IF
# endif
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND) THEN
        outwind=outTotal
      ELSE
        Vtotal=outTotal(:,1)
        CALL SCATTER_ONED_ARRAY(Vtotal, Vlocal)
        outwind(:,1)=Vlocal
        !
        Vtotal=outTotal(:,2)
        CALL SCATTER_ONED_ARRAY(Vtotal, Vlocal)
        outwind(:,2)=Vlocal
      END IF
# else
      outwind=outTotal
# endif
!     Print *, 'End of READ_GRIB_WIND'
      END SUBROUTINE
#endif
