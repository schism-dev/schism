#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_TFN(TFN)
      USE DATAPOOL
      IMPLICIT NONE
      integer, intent(out) :: TFN(2)
      CHARACTER(LEN=15) :: eTimeStr
      integer year, month, day, hour, min, sec
      REAL(rkind) eTimeDay
      eTimeDay=MAIN%TMJD
      CALL MJD2CT(eTimeDay,eTimeStr)
      CALL DATE_ConvertString2six(year, month, day, hour, min, sec, eTimeStr)
      TFN(1)=year*10000 + month*100 + day
      TFN(2)=hour*10000 + min*100   + sec
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXPORT_WIND_WW3_FORMAT
      USE DATAPOOL
      IMPLICIT NONE
      CHARACTER(LEN=13) :: TSSTR, IDSTR = 'WAVEWATCH III'
      INTEGER, PARAMETER :: UNGTYPE = 3
      LOGICAL, SAVE :: IsFirst = .TRUE.
      INTEGER :: GTYPE, NX, NY, IX, IY
      INTEGER :: FILLER(3)
      INTEGER :: TIDEFLAG = 0
      INTEGER :: TFN(2)
      INTEGER :: IP
      INTEGER TheOut
      real, allocatable :: Uwr(:,:), Vwr(:,:)
      real, allocatable :: Uwr_glob(:,:), Vwr_glob(:,:)
      CHARACTER(LEN=3)  :: IDFLD ='WND'
      FILLER(:)=0
      NX=np_total
      NY=1
      GTYPE=UNGTYPE
      CALL COMPUTE_TFN(TFN)
      allocate(Vwr(NX,NY), Uwr(NX,NY), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_export_ww3, current part')
#if defined MPI_PARALL_GRID
      DO IP=1,MNP
         Uwr(iplg(IP),1)=MySNGL(WINDXY(IP,1))
         Vwr(iplg(IP),1)=MySNGL(WINDXY(IP,2))
      END DO
      allocate(Vwr_glob(NX,NY), Uwr_glob(NX,NY), stat=istat)
      call mpi_reduce(Uwr,Uwr_glob,NP_GLOBAL,MPI_REAL, MPI_SUM,0,comm,ierr)
      call mpi_reduce(Vwr,Vwr_glob,NP_GLOBAL,MPI_REAL, MPI_SUM,0,comm,ierr)
      IF (myrank == 0) THEN
        DO IP=1,np_total
          Uwr(IP,1)=Uwr_glob(IP,1)*nwild_gb(IP)
          Vwr(IP,1)=Vwr_glob(IP,1)*nwild_gb(IP)
        END DO
      END IF
      deallocate(Vwr_glob, Uwr_glob)
#else
      DO IP=1,np_total
         Uwr(IP,1)=MySNGL(WINDXY(IP,1))
         Vwr(IP,1)=MySNGL(WINDXY(IP,2))
      END DO
#endif
#if defined MPI_PARALL_GRID
      IF (myrank == 0) THEN
#endif
        TheOut=FHNDL_EXPORT_CURR_WW3
        IF (IsFirst) THEN
          OPEN(TheOut, FILE='wind.ww3', FORM='UNFORMATTED', status='replace', action='write')
          WRITE (TheOut) IDSTR, IDFLD, NX, NY, GTYPE, FILLER(1:2), TIDEFLAG
        ELSE
          OPEN(TheOut, FILE='wind.ww3', FORM='UNFORMATTED', status='old', position='append', action='write')
        END IF
        IsFirst=.FALSE.
        WRITE (TheOut) TFN
        WRITE (TheOut) ((Vwr(IX,IY),IX=1,NX),IY=1,NY)
        WRITE (TheOut) ((Uwr(IX,IY),IX=1,NX),IY=1,NY)
        CLOSE(TheOut)
#if defined MPI_PARALL_GRID
      END IF
#endif
      deallocate(Vwr, Uwr)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXPORT_CURR_WW3_FORMAT
      USE DATAPOOL
      IMPLICIT NONE
      CHARACTER(LEN=13) :: TSSTR, IDSTR = 'WAVEWATCH III'
      INTEGER, PARAMETER :: UNGTYPE = 3
      LOGICAL, SAVE :: IsFirst = .TRUE.
      INTEGER :: GTYPE, NX, NY, IX, IY
      INTEGER :: FILLER(3)
      INTEGER :: TIDEFLAG = 0
      INTEGER :: TFN(2)
      INTEGER :: IP
      INTEGER TheOut
      real, allocatable :: Uwr(:,:), Vwr(:,:)
      real, allocatable :: Uwr_glob(:,:), Vwr_glob(:,:)
      CHARACTER(LEN=3)  :: IDFLD ='CUR'
      FILLER(:)=0
      NX=np_total
      NY=1
      GTYPE=UNGTYPE
      CALL COMPUTE_TFN(TFN)
      allocate(Vwr(NX,NY), Uwr(NX,NY), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_export_ww3, current part')
#if defined MPI_PARALL_GRID
      DO IP=1,MNP
         Uwr(iplg(IP),1)=MySNGL(CURTXY(IP,1))
         Vwr(iplg(IP),1)=MySNGL(CURTXY(IP,2))
      END DO
      allocate(Vwr_glob(NX,NY), Uwr_glob(NX,NY), stat=istat)
      call mpi_reduce(Uwr,Uwr_glob,NP_GLOBAL,MPI_REAL, MPI_SUM,0,comm,ierr)
      call mpi_reduce(Vwr,Vwr_glob,NP_GLOBAL,MPI_REAL, MPI_SUM,0,comm,ierr)
      IF (myrank == 0) THEN
        DO IP=1,np_total
          Uwr(IP,1)=Uwr_glob(IP,1)*nwild_gb(IP)
          Vwr(IP,1)=Vwr_glob(IP,1)*nwild_gb(IP)
        END DO
      END IF
      deallocate(Vwr_glob, Uwr_glob)
#else
      DO IP=1,np_total
         Uwr(IP,1)=MySNGL(CURTXY(IP,1))
         Vwr(IP,1)=MySNGL(CURTXY(IP,2))
      END DO
#endif
#if defined MPI_PARALL_GRID
      IF (myrank == 0) THEN
#endif
        TheOut=FHNDL_EXPORT_CURR_WW3
        IF (IsFirst) THEN
          OPEN(TheOut, FILE='current.ww3', FORM='UNFORMATTED', status='replace', action='write')
          WRITE (TheOut) IDSTR, IDFLD, NX, NY, GTYPE, FILLER(1:2), TIDEFLAG
        ELSE
          OPEN(TheOut, FILE='current.ww3', FORM='UNFORMATTED', status='old', position='append', action='write')
        END IF
        IsFirst=.FALSE.
        WRITE (TheOut) TFN
        WRITE (TheOut) ((Vwr(IX,IY),IX=1,NX),IY=1,NY)
        WRITE (TheOut) ((Uwr(IX,IY),IX=1,NX),IY=1,NY)
        CLOSE(TheOut)
#if defined MPI_PARALL_GRID
      END IF
#endif
      deallocate(Vwr, Uwr)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXPORT_WALV_WW3_FORMAT
      USE DATAPOOL
      IMPLICIT NONE
      CHARACTER(LEN=13) :: TSSTR, IDSTR = 'WAVEWATCH III'
      INTEGER, PARAMETER :: UNGTYPE = 3
      LOGICAL, SAVE :: IsFirst = .TRUE.
      INTEGER :: GTYPE, NX, NY, IX, IY
      INTEGER :: FILLER(3)
      INTEGER :: TIDEFLAG = 0
      INTEGER :: TFN(2)
      INTEGER :: IP
      INTEGER TheOut
      real, allocatable :: LEVwr(:,:), LEVwr_glob(:,:)
      CHARACTER(LEN=3)  :: IDFLD ='LEV'
      FILLER(:)=0
      NX=np_total
      NY=1
      GTYPE=UNGTYPE
      CALL COMPUTE_TFN(TFN)
      allocate(LEVwr(NX,NY), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_export_ww3, WALV part 1')
#if defined MPI_PARALL_GRID
      DO IP=1,MNP
         LEVwr(iplg(IP),1)=MySNGL(WATLEV(IP))
      END DO
      allocate(LEVwr_glob(NX,NY), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_export_ww3, WALV part 2')
      call mpi_reduce(LEVwr,LEVwr_glob,NP_GLOBAL,MPI_REAL, MPI_SUM,0,comm,ierr)
      IF (myrank == 0) THEN
        DO IP=1,np_total
          LEVwr(IP,1)=LEVwr_glob(IP,1)*nwild_gb(IP)
        END DO
      END IF
      deallocate(LEVwr_glob)
#else
      DO IP=1,np_total
         LEVwr(IP,1)=MySNGL(WATLEV(IP))
      END DO
#endif
#if defined MPI_PARALL_GRID
      IF (myrank == 0) THEN
#endif
        TheOut=FHNDL_EXPORT_WALV_WW3
        IF (IsFirst) THEN
          OPEN(TheOut, FILE='level.ww3', FORM='UNFORMATTED', status='replace', action='write')
          WRITE (TheOut) IDSTR, IDFLD, NX, NY, GTYPE, FILLER(1:2), TIDEFLAG
        ELSE
          OPEN(TheOut, FILE='level.ww3', FORM='UNFORMATTED', status='old', position='append', action='write')
        END IF
        IsFirst=.FALSE.
        WRITE (TheOut) TFN
        WRITE (TheOut) ((LEVwr(IX,IY),IX=1,NX),IY=1,NY)
        CLOSE(TheOut)
#if defined MPI_PARALL_GRID
      END IF
#endif
      deallocate(LEVwr)
      END SUBROUTINE
