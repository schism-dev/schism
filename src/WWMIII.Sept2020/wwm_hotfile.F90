#include "wwm_functions.h"
MODULE wwm_hotfile_mod
      USE DATAPOOL
      TYPE Subset
         integer eRankProc
         integer nbNeedEntries
         integer NPLOC
         integer, dimension(:), pointer :: ListNeedIndexFile
         integer, dimension(:), pointer :: ListNeedIndexMemory
      END TYPE Subset
      TYPE ReconstructInfo
         LOGICAL IsEasy
         integer nbNeedProc
         type(Subset), dimension(:), pointer :: ListSubset
      END TYPE ReconstructInfo
      INTEGER, parameter :: nbOned = 37
      REAL(rkind), allocatable :: ACreturn(:,:,:)
      REAL(rkind), allocatable :: VAR_ONEDreturn(:,:)
      REAL(rkind), allocatable :: VAR_ONED(:,:)
      LOGICAL :: HOTFILE_INIT = .FALSE.
      CONTAINS
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PRE_CREATE_LOCAL_HOTNAME(FILEHOT, FILERET,            &
     &  MULTIPLE, HOTSTYLE, eRank)
      IMPLICIT NONE
      character(len=*), intent(in) :: FILEHOT
      character(len=*), intent(out) :: FILERET
      integer, intent(in) :: MULTIPLE, HOTSTYLE, eRank
!
      character(len=1), parameter :: ePoint = '.'
      integer, parameter :: powerproc = 6
      character(len=6) :: eStrProc
      integer :: LPOS
      integer POSITION_BEFORE_POINT
      IF (MULTIPLE.eq.0) THEN
        LPOS=POSITION_BEFORE_POINT(FILEHOT)
        IF (HOTSTYLE == 1) THEN
          WRITE (FILERET,10) TRIM(FILEHOT(1:LPOS))
  10      FORMAT (a,'.dat')
        ELSE
          WRITE (FILERET,20) TRIM(FILEHOT(1:LPOS))
  20      FORMAT (a,'.nc')
        ENDIF
      ELSE
        LPOS=POSITION_BEFORE_POINT(FILEHOT)
        CALL GETSTRING(powerproc, eRank, eStrProc)
        IF (HOTSTYLE == 1) THEN
          WRITE (FILERET,40) TRIM(FILEHOT(1:LPOS)),eStrProc
  40      FORMAT (a,'_',a,'.dat')
        ELSE
          WRITE (FILERET,50) TRIM(FILEHOT(1:LPOS)),eStrProc
  50      FORMAT (a,'_',a,'.nc')
        ENDIF
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CREATE_LOCAL_HOTNAME(FILEHOT, FILERET, MULTIPLE, HOTSTYLE)
      IMPLICIT NONE
      character(LEN=*), intent(in) :: FILEHOT
      character(LEN=140), intent(OUT) :: FILERET
      integer, intent(in) :: MULTIPLE, HOTSTYLE
#ifdef MPI_PARALL_GRID
      CALL PRE_CREATE_LOCAL_HOTNAME(FILEHOT, FILERET,                  &
     &  MULTIPLE, HOTSTYLE, myrank+1)
#else
      FILERET=FILEHOT
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_AC_SIMPLE(FILEHOT, NPCALL, ACread, VAR_ONEDread)
      IMPLICIT NONE
      character(len=*), intent(in) :: FILEHOT
      integer, intent(in) :: NPCALL
      real(rkind), intent(inout) :: ACread(MSC, MDC, NPCALL)
      real(rkind), intent(inout) :: VAR_ONEDread(nbOned, NPCALL)
      integer :: NPLOC
      integer istat
      integer, allocatable :: IPLGloc(:)
      INTEGER :: HMNP, HMNE
      INTEGER :: HMSC, HMDC
      REAL(rkind) :: HFRLOW, HFRHIGH
      OPEN(HOTIN%FHNDL, FILE = TRIM(FILEHOT), STATUS = 'OLD', FORM = 'UNFORMATTED')
      READ(HOTIN%FHNDL) HMNP, HMNE
      READ(HOTIN%FHNDL) HMSC, HMDC, HFRLOW, HFRHIGH
      IF ( HMNP .NE. NP_TOTAL .OR. HMNE .NE. NE_TOTAL .OR.          &
     &     HMSC .NE. MSC      .OR. HFRLOW .NE. FRLOW .OR.           &
     &    HFRHIGH .NE. FRHIGH ) THEN
        CALL WWM_ABORT('THE HOTFILE GEOMETRY DOES NOT FIT THE INPUT FILE')
      ENDIF
      READ(HOTIN%FHNDL) NPLOC
      allocate(IPLGloc(NPLOC), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 1')
      READ(HOTIN%FHNDL) IPLGloc
      deallocate(IPLGloc)
!todo ordering      
      READ(HOTIN%FHNDL) ACread
      READ(HOTIN%FHNDL) VAR_ONEDread
      CLOSE(HOTIN%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DETERMINE_NUMBER_PROC(FILEHOT, HOTSTYLE, nbProc)
      IMPLICIT NONE
      character(len=*), intent(in) :: FILEHOT
      integer, intent(in) :: HOTSTYLE
      integer, intent(out) :: nbProc
      INTEGER :: MULTIPLE = 1
      character(len=140) :: FILERET
      integer :: iRankTest
      logical :: test
      iRankTest=1
      DO
        CALL PRE_CREATE_LOCAL_HOTNAME(FILEHOT, FILERET, MULTIPLE, HOTSTYLE, iRankTest)
        INQUIRE(FILE=TRIM(FILERET), EXIST=test)
        IF (test.eqv..false.) THEN
          nbProc=iRankTest-1
          EXIT
        END IF
        iRankTest=iRankTest+1
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_IPLG(HOTSTYLE, FILEHOT, eRank, NPLOC, IPLG)
#ifdef NCDF
      USE NETCDF
#endif
      IMPLICIT NONE
      character (len = *), parameter :: CallFct="READ_IPLG"
      character(len=*), intent(in) :: FILEHOT
      integer, intent(IN) :: HOTSTYLE, eRank
      integer, intent(out) :: NPLOC, IPLG(np_total)
#ifdef NCDF
      integer :: iret, ncid, mnp_dims
#endif
      INTEGER :: HMNP, HMNE
      INTEGER :: HMSC, HMDC
      INTEGER, allocatable :: IPLGin(:)
#ifdef NCDF
      INTEGER :: iplg_id
#endif
      REAL(rkind) :: HFRLOW, HFRHIGH
      integer :: MULTIPLE, istat
      character(len=140) :: FILERET
      MULTIPLE=1
      CALL PRE_CREATE_LOCAL_HOTNAME(FILEHOT, FILERET, MULTIPLE, HOTSTYLE, eRank)
      IF (HOTSTYLE.eq.1) THEN
        OPEN(HOTIN%FHNDL, FILE = TRIM(FILERET), STATUS = 'OLD', FORM = 'UNFORMATTED')
        READ(HOTIN%FHNDL) HMNP, HMNE
        READ(HOTIN%FHNDL) HMSC, HMDC, HFRLOW, HFRHIGH
        IF ( HMNP .NE. NP_TOTAL .OR. HMNE .NE. NE_TOTAL .OR.            &
     &       HMSC .NE. MSC      .OR. HFRLOW .NE. FRLOW .OR.             &
     &       HFRHIGH .NE. FRHIGH ) THEN
          CALL WWM_ABORT('THE HOTFILE GEOMETRY DOES NOT FIT THE INPUT FILE')
        ENDIF
        READ(HOTIN%FHNDL) NPLOC
        allocate(IPLGin(NPLOC), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 2')
        READ(HOTIN%FHNDL) IPLGin
        IPLG(1:NPLOC)=IPLGin
        deallocate(IPLGin)
        CLOSE(HOTIN%FHNDL)
      ELSE
#ifdef NCDF
        iret=nf90_open(TRIM(FILERET), nf90_nowrite, ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
        iret=nf90_inq_dimid(ncid,"mnp",mnp_dims)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)
        iret=nf90_inquire_dimension(ncid, mnp_dims, len=nploc)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)
        allocate(IPLGin(NPLOC), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 3')
        iret=nf90_inq_varid(ncid, "iplg", iplg_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
        iret=nf90_get_var(ncid, iplg_id, IPLGin, start=(/1/), count=(/NPLOC/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
        iret=nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
        IPLG(1:NPLOC)=IPLGin
        deallocate(IPLGin)
#endif
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEALLOCATE_ReconsArr(eRecons)
      IMPLICIT NONE
      type(ReconstructInfo), intent(inout) :: eRecons
      integer :: nbNeedProc, iProc
      IF (eRecons % IsEasy.eqv..false.) THEN
        nbNeedProc=eRecons % nbNeedProc
        DO iProc=1,nbNeedProc
          DEALLOCATE(eRecons % ListSubset(iProc) % ListNeedIndexFile)
          DEALLOCATE(eRecons % ListSubset(iProc) % ListNeedIndexMemory)
        END DO
        DEALLOCATE(eRecons % ListSubset)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DETERMINE_NEEDED_HOTFILES(HOTSTYLE, FILEHOT, eRecons)
      IMPLICIT NONE
      INTEGER, intent(in) :: HOTSTYLE
      character(len=*), intent(in) :: FILEHOT
      type(ReconstructInfo), intent(inout) :: eRecons
      character(len=140) :: errmsg
      integer, allocatable :: IPLGtot(:)
      integer, allocatable :: eStatus(:)
      integer, allocatable :: ListAttained(:)
      integer :: eDiff, I, iProc, idx, eIdx, nbProc, IP
      integer :: nbNeedProc, nbF, NPLOC, nbZero, idxB
      integer istat
#ifndef MPI_PARALL_GRID
      integer, allocatable :: iplg(:)
      integer :: nproc, myrank
      nproc=1
      myrank=0
      allocate(iplg(MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 4')
      DO IP=1,MNP
        iplg(IP)=IP
      END DO
#endif
      allocate(IPLGtot(np_total), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 5')
      CALL DETERMINE_NUMBER_PROC(FILEHOT, HOTSTYLE, nbProc)
      IF (nbProc.eq.nproc) THEN
        CALL READ_IPLG(HOTSTYLE, FILEHOT, myrank+1, NPLOC, IPLGtot)
        IF (NPLOC.eq.MNP) THEN
          eDiff=0
          DO IP=1,MNP
            eDiff=eDiff + abs(IPLG(IP) - IPLGtot(IP))
          END DO
          IF (eDiff.eq.0) THEN
            eRecons % IsEasy = .TRUE.
            deallocate(IPLGtot)
            RETURN
          END IF
        END IF
      END IF
      eRecons % IsEasy = .FALSE.
      allocate(eStatus(np_total), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 6')
      eStatus=0
      DO IP=1,MNP
        eStatus(IPLG(IP))=IP
      END DO
      nbNeedProc=0
      IF (myrank.eq.0) THEN
        allocate(ListAttained(np_total), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 7')
        ListAttained=0
      ENDIF
      DO iProc=1,nbProc
        CALL READ_IPLG(HOTSTYLE, FILEHOT, iProc, NPLOC, IPLGtot)
        nbF=0
        DO I=1,NPLOC
          IF (eStatus(IPLGtot(I)).gt.0) THEN
            nbF=nbF+1
          END IF
          IF (myrank.eq.0) THEN
            ListAttained(IPLGtot(I))=1
          ENDIF
        END DO
        IF (nbF.gt.0) THEN
          nbNeedProc=nbNeedProc+1
        END IF
      END DO
      IF (myrank.eq.0) THEN
        nbZero=0
        DO IP=1,np_total
          IF (ListAttained(IP).eq.0) THEN
            nbZero=nbZero+1
          ENDIF
        END DO
        IF (nbZero.gt.0) THEN
          WRITE(errmsg, *) 'Not enough data nbProc=', nbProc, ' nbMissedMode=', nbZero
          CALL WWM_ABORT(errmsg)
        ENDIF
        DEALLOCATE(ListAttained)
      ENDIF
      eRecons % nbNeedProc=nbNeedProc
      allocate(eRecons % ListSubset(nbNeedProc), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 8')
      idx=0
      DO iProc=1,nbProc
        CALL READ_IPLG(HOTSTYLE, FILEHOT, iProc, NPLOC, IPLGtot)
        nbF=0
        DO IP=1,NPLOC
          IF (eStatus(IPLGtot(IP)).gt.0) THEN
            nbF=nbF+1
          END IF
        END DO
        IF (nbF.gt.0) THEN
          idx=idx+1
          eRecons % ListSubset(idx) % eRankProc=iProc
          eRecons % ListSubset(idx) % nbNeedEntries=nbF
          eRecons % ListSubset(idx) % NPLOC=NPLOC
          ALLOCATE(eRecons % ListSubset(idx) % ListNeedIndexFile(nbF), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 9')
          ALLOCATE(eRecons % ListSubset(idx) % ListNeedIndexMemory(nbF), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 10')
          idxB=0
          DO I=1,NPLOC
            eIdx=IPLGtot(I)
            IP=eStatus(eIdx)
            IF (IP.gt.0) THEN
              idxB=idxB+1
              eRecons % ListSubset(idx) % ListNeedIndexFile(idxB)=I
              eRecons % ListSubset(idx) % ListNeedIndexMemory(idxB)=IP
            END IF
          END DO
        END IF
      END DO
      deallocate(eStatus, IPLGtot)
#ifndef MPI_PARALL_GRID
      deallocate(iplg)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INPUT_HOTFILE
      IMPLICIT NONE
      integer IP
      allocate(VAR_ONED(nbOned, MNP), stat=istat)
      IF (HOTSTYLE_IN == 1) THEN
        CALL INPUT_HOTFILE_BINARY
#ifdef NCDF
      ELSE IF (HOTSTYLE_IN == 2) THEN
        CALL INPUT_HOTFILE_NETCDF
#endif
      ELSE
        CALL WWM_ABORT('Wrong choice of HOTSTYLE_IN')
      END IF
      DO IP=1,MNP
        WINDXY(IP,1)=VAR_ONED(1 ,IP)
        WINDXY(IP,2)=VAR_ONED(2 ,IP)
        PRESSURE(IP)=VAR_ONED(3 ,IP)
        DVWIND(IP,1)=VAR_ONED(4 ,IP)
        DVWIND(IP,2)=VAR_ONED(5 ,IP)
        CURTXY(IP,1)=VAR_ONED(6 ,IP)
        CURTXY(IP,2)=VAR_ONED(7 ,IP)
        DVCURT(IP,1)=VAR_ONED(8 ,IP)
        DVCURT(IP,2)=VAR_ONED(9 ,IP)
        DDEP(IP,1)=VAR_ONED(10,IP)
        DDEP(IP,2)=VAR_ONED(11,IP)
        DCUX(IP,1)=VAR_ONED(12,IP)
        DCUX(IP,2)=VAR_ONED(13,IP)
        DCUY(IP,1)=VAR_ONED(14,IP)
        DCUY(IP,2)=VAR_ONED(15,IP)
        WATLEV(IP)=VAR_ONED(16,IP)
        WATLEVOLD(IP)=VAR_ONED(17,IP)
        WLDEP(IP)=VAR_ONED(18,IP)
        DEPDT(IP)=VAR_ONED(19,IP)
        QBLOCAL(IP)=VAR_ONED(20,IP)
        DISSIPATION(IP)=VAR_ONED(21,IP)
        AIRMOMENTUM(IP)=VAR_ONED(22,IP)
        UFRIC(IP)=VAR_ONED(23,IP)
        ALPHA_CH(IP)=VAR_ONED(24,IP)
        TAUW(IP)=VAR_ONED(25,IP)
        TAUTOT(IP)=VAR_ONED(26,IP)
        TAUWX(IP)=VAR_ONED(27,IP)
        TAUWY(IP)=VAR_ONED(28,IP)
        TAUHF(IP)=VAR_ONED(29,IP)
        Z0(IP)=VAR_ONED(30,IP)
        CD(IP)=VAR_ONED(31,IP)
        USTDIR(IP)=VAR_ONED(32,IP)
        RSXX(IP)=VAR_ONED(33,IP)
        RSXY(IP)=VAR_ONED(34,IP)
        RSYY(IP)=VAR_ONED(35,IP)
        FORCEXY(IP,1)=VAR_ONED(36,IP)
        FORCEXY(IP,2)=VAR_ONED(37,IP)
      END DO
      deallocate(VAR_ONED)
      END SUBROUTINE
!**********************************************************************
!*  Now the other side of the story. The output.                      *
!*  Everything in reverse                                             *
!**********************************************************************
      SUBROUTINE INIT_HOTFILE_OUTPUT
      IMPLICIT NONE
#if defined NCDF && defined MPI_PARALL_GRID
      integer :: ListFirst(nproc)
      integer MNPloc, iProc, IP, IP_glob
      integer, allocatable :: dspl_ac(:), dspl_var_oned(:)
      ListFirst=0
      DO iProc=2,nproc
        ListFirst(iProc)=ListFirst(iProc-1) + ListMNP(iProc-1)
      END DO
      IF (myrank .eq. 0) THEN
        allocate(ac2_hot_rqst(nproc-1), ac2_hot_stat(MPI_STATUS_SIZE,nproc-1), ac2_hot_type(nproc-1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('INIT_HOTFILE_OUTPUT, allocate error 1')
        allocate(var_oned_hot_rqst(nproc-1), var_oned_hot_stat(MPI_STATUS_SIZE,nproc-1), var_oned_hot_type(nproc-1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('INIT_HOTFILE_OUTPUT, allocate error 2')
        DO iProc=2,nproc
          MNPloc=ListMNP(iProc)
          allocate(dspl_ac(MNPloc), dspl_var_oned(MNPloc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('error in IOBPtotal allocate')
          DO IP=1,MNPloc
            IP_glob=ListIPLG(IP+ListFirst(iProc))
            dspl_ac(IP)=MSC*MDC*(IP_glob-1)
            dspl_var_oned(IP)=nbOned*(IP_glob-1)
          END DO
          call mpi_type_create_indexed_block(MNPloc,MSC*MDC,dspl_ac,rtype,ac2_hot_type(iProc-1), ierr)
          call mpi_type_commit(ac2_hot_type(iProc-1), ierr)
          call mpi_type_create_indexed_block(MNPloc,nbOned,dspl_var_oned,rtype,var_oned_hot_type(iProc-1), ierr)
          call mpi_type_commit(var_oned_hot_type(iProc-1), ierr)
          deallocate(dspl_ac, dspl_var_oned)
        END DO
      END IF
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SETUP_RETURN_AC_VARONED
      IMPLICIT NONE
#if defined NCDF && defined MPI_PARALL_GRID
      integer iProc, IP, IPglob
      IF (myrank .eq. 0) THEN
        WRITE(STAT%FHNDL, *) 'Before allocation of ACreturn'
        allocate(ACreturn(MSC,MDC,np_global), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 16')
        DO iProc=2,nproc
          call mpi_irecv(ACreturn,1,ac2_hot_type(iProc-1),iProc-1,8123,comm,ac2_hot_rqst(iProc-1),ierr)
        END DO
        DO IP=1,NP_RES
          IPglob=iplg(IP)
          ACreturn(:,:,IPglob)=AC2(:,:,IP)
        END DO
        IF (nproc > 1) THEN
          call mpi_waitall(nproc-1, ac2_hot_rqst, ac2_hot_stat,ierr)
        END IF
      ELSE
        CALL MPI_SEND(AC2, MSC*MDC*NP_RES, rtype, 0, 8123, comm, ierr)
      END IF
      IF (myrank .eq. 0) THEN
        WRITE(STAT%FHNDL, *) 'Before allocation of VAR_ONEDreturnn'
        allocate(VAR_ONEDreturn(nbOned,np_global), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 16')
        DO iProc=2,nproc
          call mpi_irecv(VAR_ONEDreturn,1,var_oned_hot_type(iProc-1),iProc-1,8124,comm,var_oned_hot_rqst(iProc-1),ierr)
        END DO
        DO IP=1,NP_RES
          IPglob=iplg(IP)
          VAR_ONEDreturn(:,IPglob)=VAR_ONED(:,IP)
        END DO
        IF (nproc > 1) THEN
          call mpi_waitall(nproc-1,var_oned_hot_rqst,var_oned_hot_stat,ierr)
        END IF
      ELSE
        CALL MPI_SEND(VAR_ONED, nbOned*NP_RES, rtype, 0, 8124, comm, ierr)
      END IF
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OUTPUT_HOTFILE
      IMPLICIT NONE
      integer IP
      allocate(VAR_ONED(nbOned, MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 11')
      DO IP=1,MNP
        VAR_ONED(1 ,IP)=WINDXY(IP,1)
        VAR_ONED(2 ,IP)=WINDXY(IP,2)
        VAR_ONED(3 ,IP)=PRESSURE(IP)
        VAR_ONED(4 ,IP)=DVWIND(IP,1)
        VAR_ONED(5 ,IP)=DVWIND(IP,2)
        VAR_ONED(6 ,IP)=CURTXY(IP,1)
        VAR_ONED(7 ,IP)=CURTXY(IP,2)
        VAR_ONED(8 ,IP)=DVCURT(IP,1)
        VAR_ONED(9 ,IP)=DVCURT(IP,2)
        VAR_ONED(10,IP)=DDEP(IP,1)
        VAR_ONED(11,IP)=DDEP(IP,2)
        VAR_ONED(12,IP)=DCUX(IP,1)
        VAR_ONED(13,IP)=DCUX(IP,2)
        VAR_ONED(14,IP)=DCUY(IP,1)
        VAR_ONED(15,IP)=DCUY(IP,2)
        VAR_ONED(16,IP)=WATLEV(IP)
        VAR_ONED(17,IP)=WATLEVOLD(IP)
        VAR_ONED(18,IP)=WLDEP(IP)
        VAR_ONED(19,IP)=DEPDT(IP)
        VAR_ONED(20,IP)=QBLOCAL(IP)
        VAR_ONED(21,IP)=DISSIPATION(IP)
        VAR_ONED(22,IP)=AIRMOMENTUM(IP)
        VAR_ONED(23,IP)=UFRIC(IP)
        VAR_ONED(24,IP)=ALPHA_CH(IP)
        VAR_ONED(25,IP)=TAUW(IP)
        VAR_ONED(26,IP)=TAUTOT(IP)
        VAR_ONED(27,IP)=TAUWX(IP)
        VAR_ONED(28,IP)=TAUWY(IP)
        VAR_ONED(29,IP)=TAUHF(IP)
        VAR_ONED(30,IP)=Z0(IP)
        VAR_ONED(31,IP)=CD(IP)
        VAR_ONED(32,IP)=USTDIR(IP)
        VAR_ONED(33,IP)=RSXX(IP)
        VAR_ONED(34,IP)=RSXY(IP)
        VAR_ONED(35,IP)=RSYY(IP)
        VAR_ONED(36,IP)=FORCEXY(IP,1)
        VAR_ONED(37,IP)=FORCEXY(IP,2)
      END DO
#ifdef MPI_PARALL_GRID
      IF (.NOT. HOTFILE_INIT) THEN
        HOTFILE_INIT=.TRUE.
        CALL INIT_HOTFILE_OUTPUT
      END IF
      IF (MULTIPLEOUT_HOT.eq.0) THEN
        CALL SETUP_RETURN_AC_VARONED
      END IF
#endif
      IF (HOTSTYLE_OUT == 1) THEN
        CALL OUTPUT_HOTFILE_BINARY
#ifdef NCDF
      ELSE IF (HOTSTYLE_OUT == 2) THEN
        CALL OUTPUT_HOTFILE_NETCDF
#endif
      ELSE
        CALL WWM_ABORT('Wrong choice of HOTSTYLE_OUT')
      END IF
#ifdef MPI_PARALL_GRID
      IF ((MULTIPLEOUT_HOT.eq.0).and.(myrank.eq.0)) THEN
        WRITE(STAT%FHNDL, *) 'Before deallocation of ACreturn'
        deallocate(ACreturn)
        WRITE(STAT%FHNDL, *) 'Before deallocation of VAR_ONEDreturn'
        deallocate(VAR_ONEDreturn)
      END IF
#endif
      deallocate(VAR_ONED)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INPUT_HOTFILE_BINARY
      IMPLICIT NONE
      INTEGER idxFil, idxMem, iProc, istat
      REAL(rkind), ALLOCATABLE :: ACinB(:,:,:)
      REAL(rkind), ALLOCATABLE :: VAR_ONED_B(:,:)
      character(len=140) :: FILERET
      type(ReconstructInfo) :: eRecons
      integer :: nbF, NPLOC, eRank, I
#ifdef MPI_PARALL_GRID
      integer IP
#endif
      IF (MULTIPLEIN_HOT.eq.0) THEN
#ifdef MPI_PARALL_GRID
        ALLOCATE(ACinB(MSC,MDC,NP_GLOBAL), VAR_ONED_B(nbOned, NP_GLOBAL), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 11')
        CALL READ_AC_SIMPLE(HOTIN%FNAME, NP_GLOBAL, ACinB, VAR_ONED_B)
        DO IP=1,MNP
          AC2(:,:,IP)=ACinB(:,:,iplg(IP))
          VAR_ONED(:,IP)=VAR_ONED_B(:,iplg(IP))
        ENDDO
        DEALLOCATE(ACinB, VAR_ONED_B)
#else
        CALL READ_AC_SIMPLE(HOTIN%FNAME, MNP, AC2, VAR_ONED)
#endif
      ELSE
        CALL DETERMINE_NEEDED_HOTFILES(HOTSTYLE_IN, TRIM(HOTIN%FNAME), eRecons)
        IF (eRecons % IsEasy) THEN
          CALL CREATE_LOCAL_HOTNAME(HOTIN%FNAME, FILERET, MULTIPLEIN_HOT, HOTSTYLE_IN)
          CALL READ_AC_SIMPLE(FILERET, MNP, AC2, VAR_ONED)
        ELSE
          DO iProc=1,eRecons % nbNeedProc
            eRank=eRecons % ListSubset(iProc) % eRankProc
            nbF=eRecons % ListSubset(iProc) % nbNeedEntries
            NPLOC=eRecons % ListSubset(iProc) % NPLOC
            CALL PRE_CREATE_LOCAL_HOTNAME(HOTIN%FNAME, FILERET, MULTIPLEIN_HOT, HOTSTYLE_IN, eRank)
            allocate(ACinB(MSC, MDC,NPLOC), VAR_ONED_B(nbOned, NPLOC), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 12')
            CALL READ_AC_SIMPLE(FILERET, NPLOC, ACinB, VAR_ONED_B)
            DO I=1,nbF
              idxFil=eRecons % ListSubset(iProc) % ListNeedIndexFile(I)
              idxMem=eRecons % ListSubset(iProc) % ListNeedIndexMemory(I)
              AC2(:,:,idxMem)=ACinB(:,:,idxFil)
              VAR_ONED(:,idxMem)=VAR_ONED_B(:,idxFil)
            END DO
            deallocate(ACinB, VAR_ONED_B)
          END DO
        END IF
        CALL DEALLOCATE_ReconsArr(eRecons)
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OUTPUT_HOTFILE_BINARY
      IMPLICIT NONE
#ifdef MPI_PARALL_GRID
      include 'mpif.h'
#endif
      CHARACTER(len=140) :: FILERET
#ifdef MPI_PARALL_GRID
      integer IP, IS, ID, istat
      REAL(rkind), allocatable :: VALB(:), VALB_SUM(:)
#endif
      CALL CREATE_LOCAL_HOTNAME(HOTOUT%FNAME, FILERET, MULTIPLEOUT_HOT, HOTSTYLE_OUT)
      OPEN(HOTOUT%FHNDL, FILE = TRIM(FILERET), STATUS = 'UNKNOWN',  FORM = 'UNFORMATTED')
      WRITE(HOTOUT%FHNDL) NP_TOTAL, NE_TOTAL
      WRITE(HOTOUT%FHNDL) MSC, MDC, FRLOW, FRHIGH
#ifndef MPI_PARALL_GRID
      WRITE(HOTOUT%FHNDL) AC2
      WRITE(HOTOUT%FHNDL) VAR_ONED
#else
      IF (MULTIPLEOUT_HOT.eq.0) THEN
        WRITE(HOTOUT%FHNDL) ACreturn
        WRITE(HOTOUT%FHNDL) VAR_ONEDreturn
      ELSE
        WRITE(HOTOUT%FHNDL) MNP
        WRITE(HOTOUT%FHNDL) IPLG
        WRITE(HOTOUT%FHNDL) AC2
        WRITE(HOTOUT%FHNDL) VAR_ONED
      ENDIF
#endif
      CLOSE(HOTOUT%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef NCDF
      SUBROUTINE INPUT_HOTFILE_NETCDF
      USE NETCDF
      IMPLICIT NONE
# ifdef MPI_PARALL_GRID
      INTEGER :: IP, ID
      REAL(rkind) :: ACLOC(MSC,MDC)
      REAL(rkind) :: VARLOC(nbOned)
# endif
      INTEGER :: NPLOC, eRank, I
      character (len = *), parameter :: CallFct="INPUT_HOTFILE_NETCDF"
      REAL(rkind), ALLOCATABLE :: ACinB(:,:,:)
      REAL(rkind), ALLOCATABLE :: VAR_ONED_B(:,:)
      type(ReconstructInfo) :: eRecons
      character(len=140) :: FILERET
      INTEGER :: iret, ncid, ac_id, var_oned_id
      INTEGER :: nbF, iProc, idxFil, idxMem
      integer istat
      IF (MULTIPLEIN_HOT.eq.0) THEN
# ifdef MPI_PARALL_GRID
        iret=nf90_open(TRIM(HOTIN%FNAME), nf90_nowrite, ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
        iret=nf90_inq_varid(ncid, "ac", ac_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)
        iret=nf90_inq_varid(ncid, "var_oned", var_oned_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)
        DO IP=1,MNP
          iret=nf90_get_var(ncid,ac_id,ACLOC, start=(/1,1,iplg(IP),IHOTPOS_IN/), count = (/MSC, MDC, 1, 1 /))
          IF (iret /= 0) THEN
            Print *, 'This time send direcly your bug to Mathieu.Dutour@gmail.com'
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
          END IF
          AC2(:,:,IP)=ACLOC
          iret=nf90_get_var(ncid,var_oned_id,VARLOC, start=(/1,iplg(IP),IHOTPOS_IN/), count = (/nbOned, 1, 1 /))
          IF (iret /= 0) THEN
            Print *, 'Same story. Send direcly your bug to Mathieu.Dutour@gmail.com'
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
          END IF
          VAR_ONED(:,IP)=VARLOC
        END DO
        iret=nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
# else
        iret=nf90_open(TRIM(HOTIN%FNAME), nf90_nowrite, ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
        iret=nf90_inq_varid(ncid, "ac", ac_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, iret)
        iret=nf90_inq_varid(ncid, "var_oned", var_oned_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, iret)
        iret=nf90_get_var(ncid,ac_id,AC2, start=(/1,1,1,IHOTPOS_IN/), count=(/MSC,MDC,MNP,1/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, iret)
        iret=nf90_get_var(ncid,var_oned_id,VAR_ONED, start=(/1,1,IHOTPOS_IN/), count=(/nbOned, MNP, 1/))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, iret)
        iret=nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, iret)
# endif
      ELSE
        CALL DETERMINE_NEEDED_HOTFILES(HOTSTYLE_IN, TRIM(HOTIN%FNAME), eRecons)
        IF (eRecons % IsEasy) THEN
          CALL CREATE_LOCAL_HOTNAME(HOTIN%FNAME, FILERET, MULTIPLEIN_HOT, HOTSTYLE_IN)
          iret=nf90_open(FILERET, nf90_nowrite, ncid)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, iret)
          iret=nf90_inq_varid(ncid, "ac", ac_id)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, iret)
          iret=nf90_inq_varid(ncid, "var_oned", var_oned_id)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, iret)
          iret=nf90_get_var(ncid,ac_id,AC2, start=(/1,1,1,IHOTPOS_IN/),  count = (/MSC, MDC, MNP, 1 /))
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 15, iret)
          iret=nf90_get_var(ncid,var_oned_id,VAR_ONED, start=(/1,1,IHOTPOS_IN/), count=(/nbOned, MNP, 1/))
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 16, iret)
          iret=nf90_close(ncid)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 17, iret)
        ELSE
          DO iProc=1,eRecons % nbNeedProc
            eRank=eRecons % ListSubset(iProc) % eRankProc
            nbF=eRecons % ListSubset(iProc) % nbNeedEntries
            NPLOC=eRecons % ListSubset(iProc) % NPLOC
            CALL PRE_CREATE_LOCAL_HOTNAME(HOTIN%FNAME, FILERET, MULTIPLEIN_HOT, HOTSTYLE_IN, eRank)
            iret=nf90_open(TRIM(FILERET), nf90_nowrite, ncid)
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 17, iret)
            allocate(ACinB(MSC,MDC,NPLOC), VAR_ONED_B(nbOned, NPLOC), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 15')
            iret=nf90_inq_varid(ncid, "ac", ac_id)
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 18, iret)
            iret=nf90_inq_varid(ncid, "var_oned", var_oned_id)
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 19, iret)
            iret=nf90_get_var(ncid,ac_id,ACinB, start=(/1,1,1,IHOTPOS_IN/), count=(/MSC, MDC, NPLOC, 1 /))
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 20, iret)
            iret=nf90_get_var(ncid,var_oned_id,VAR_ONED_B, start=(/1,1,IHOTPOS_IN/), count=(/nbOned, NPLOC, 1 /))
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 21, iret)
            iret=nf90_close(ncid)
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 22, iret)
            DO I=1,nbF
              idxFil=eRecons % ListSubset(iProc) % ListNeedIndexFile(I)
              idxMem=eRecons % ListSubset(iProc) % ListNeedIndexMemory(I)
              AC2(:,:,idxMem)=ACinB(:,:,idxFil)
              VAR_ONED(:,idxMem)=VAR_ONED_B(:,idxFil)
            END DO
            deallocate(ACinB)
          END DO
        ENDIF
        CALL DEALLOCATE_ReconsArr(eRecons)
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_HOTFILE_PART_1(FILERET, nbTime, MULTIPLEOUT_W, GRIDWRITE_W, IOBPD_HISTORY_W, np_write, ne_write)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
      character(len=140), intent(in) :: FILERET
      integer, intent(in) :: nbTime, MULTIPLEOUT_W
      logical, intent(in) :: GRIDWRITE_W, IOBPD_HISTORY_W
      integer, intent(in) :: np_write, ne_write
      !
      character (len = *), parameter :: CallFct="WRITE_HOTFILE_PART_1"
      character (len = *), parameter :: UNITS = "units"
      integer iret, ncid
      integer nboned_dims, nfreq_dims, ndir_dims, ntime_dims, mnp_dims
      integer ac_id
      iret = nf90_create(FILERET, NF90_CLOBBER, ncid)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)

      CALL WRITE_NETCDF_HEADERS_1(ncid, nbTime, MULTIPLEOUT_W, GRIDWRITE_W, IOBPD_HISTORY_W, np_write, ne_write)

      iret=nf90_def_dim(ncid, "nboned", nbOned, nboned_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)

      iret=nf90_inq_dimid(ncid, "mnp", mnp_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)

      iret=nf90_inq_dimid(ncid, "nfreq", nfreq_dims)
      CALL GENERIC_NETCDF_ERROR_WWM_CLEAR(ncid, CallFct, 4, iret)

      iret=nf90_inq_dimid(ncid, "ndir", ndir_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)

      iret=nf90_inq_dimid(ncid, 'ocean_time', ntime_dims)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)

      iret=nf90_def_var(ncid,"ac",NF90_RUNTYPE,(/ nfreq_dims, ndir_dims, mnp_dims, ntime_dims/),ac_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, iret)

      iret=nf90_def_var(ncid,"var_oned",NF90_RUNTYPE,(/ nboned_dims, mnp_dims, ntime_dims/),ac_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, iret)

      iret=nf90_put_att(ncid,ac_id,UNITS,'unknown')
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, iret)

      iret=nf90_close(ncid)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, iret)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_HOTFILE_PART_2(FILERET, eTimeDay, POS, np_write, ACwrite, VAR_ONEDwrite)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
      character(len=140), intent(in) :: FILERET
      real(rkind), intent(in) :: eTimeDay
      integer, intent(in) :: POS, np_write
      real(rkind), intent(in) :: ACwrite(MSC,MDC,np_write), VAR_ONEDwrite(nbOned, np_write)
      character (len = *), parameter :: CallFct="WRITE_HOTFILE_PART_2"
      integer iret, ncid, var_oned_id, ac_id
      iret=nf90_open(FILERET, nf90_write, ncid)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
      CALL WRITE_NETCDF_TIME(ncid, POS, eTimeDay)
      iret=nf90_inq_varid(ncid, "ac", ac_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)
      iret=nf90_inq_varid(ncid, "var_oned", var_oned_id)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)
      iret=nf90_put_var(ncid,ac_id,ACwrite,start=(/1, 1, 1, POS/), count=(/ MSC, MDC, np_write, 1 /))
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
      iret=nf90_put_var(ncid,var_oned_id,VAR_ONEDwrite,start=(/1, 1, POS/), count=(/ nbOned, np_write, 1 /))
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
      iret=nf90_close(ncid)
      CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OUTPUT_HOTFILE_NETCDF
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
!# ifdef MPI_PARALL_GRID
!      include 'mpif.h'
!# endif
      character (len = *), parameter :: CallFct="OUTPUT_HOTFILE_NETCDF"
      INTEGER :: POS
      integer :: iret, ncid, ntime_dims, mnp_dims, nfreq_dims, ndir_dims
      integer :: ac_id, nboned_dims, var_oned_id
      integer :: nbTime
      REAL(rkind)  :: eTimeDay
      character(len=140) :: FILERET
      integer np_write, ne_write
#ifdef MPI_PARALL_GRID
      integer ID, IS, IP !, ISTAT
#endif
#ifdef MPI_PARALL_GRID
      IF (MULTIPLEOUT_HOT.eq.0) THEN
        np_write=np_global
        ne_write=ne_global
      ELSE
        np_write=MNP
        ne_write=MNE
      ENDIF
#else
      np_write=np_global
      ne_write=ne_global
#endif
      CALL CREATE_LOCAL_HOTNAME(HOTOUT%FNAME, FILERET, MULTIPLEOUT_HOT, HOTSTYLE_OUT)
      IF (IDXHOTOUT.eq.0) THEN
!$OMP MASTER
        IF (LCYCLEHOT) THEN
          nbTime=2
        ELSE
          nbTime=-1
        END IF
        IF (WriteOutputProcess_hot) THEN
          CALL WRITE_HOTFILE_PART_1(FILERET, nbTime, MULTIPLEOUT_HOT, GRIDWRITE, IOBPD_HISTORY, np_write, ne_write)
        END IF
        CALL WRITE_NETCDF_HEADERS_2(FILERET, MULTIPLEOUT_HOT, WriteOutputProcess_hot, GRIDWRITE, np_write, ne_write)
!$OMP END MASTER
      END IF
      IF (WriteOutputProcess_hot) THEN
        IF (LCYCLEHOT) THEN
          POS=mod(IDXHOTOUT,2)+1
        ELSE
          POS=IDXHOTOUT+1
        END IF
        eTimeDay=MAIN%TMJD
        !
#ifdef MPI_PARALL_GRID
        IF (MULTIPLEOUT_HOT.eq.0) THEN
          CALL WRITE_HOTFILE_PART_2(FILERET, eTimeDay, POS, np_global, ACreturn, VAR_ONEDreturn)
        ELSE
          CALL WRITE_HOTFILE_PART_2(FILERET, eTimeDay, POS, MNP, AC2, VAR_ONED)
        ENDIF
#else
        CALL WRITE_HOTFILE_PART_2(FILERET, eTimeDay, POS, np_global, AC2, VAR_ONED)
#endif
      ENDIF
      IDXHOTOUT=IDXHOTOUT+1
      END SUBROUTINE
#endif
END MODULE wwm_hotfile_mod
