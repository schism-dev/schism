!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SINGLE_READ_SPATIAL_GRID_TOTAL(eGrid, DimMode, LVAR1D, Lsphe, eGRD, iGridType)
      USE DATAPOOL, only : rkind, istat, GridInformation, FILEDEF, DBG
#ifdef NCDF
      USE NETCDF
#endif
      IMPLICIT NONE
      type(GridInformation), intent(out) :: eGrid
      integer, intent(in) :: DimMode
      logical, intent(in) :: LVAR1D
      logical, intent(in) :: Lsphe
      type(FILEDEF), intent(in) :: eGRD
      integer, intent(in) :: iGridType
      !
      INTEGER :: I, IP, IE, ITMP, JTMP
      REAL(rkind)  :: XPDTMP, YPDTMP, ZPDTMP
      REAL(rkind) DXP1, DXP2, DXP3, DYP1, DYP2, DYP3
      INTEGER KTMP, LTMP, MTMP, NTMP, OTMP
      CHARACTER(LEN=100) :: RHEADER
#ifdef NCDF
      INTEGER :: ncid, dimidsB(2), dimidsA(1)
      character (len=20) :: MNEstr, MNPstr
      INTEGER var_id1, var_id2, var_id
      character (len = *), parameter :: CallFct="SINGLE_READ_SPATIAL_GRID_TOTAL"
#endif
      CALL TEST_FILE_EXIST_DIE('Missing grid file : ', eGRD % FNAME)
      SELECT CASE (DimMode)
        CASE (1)
           !
           ! This code is broken.
           ! The np_total has not been assigned here
           !
          OPEN(eGRD % FHNDL, FILE = eGRD % FNAME, STATUS = 'OLD')
          allocate(eGrid % XPtotal(eGrid % np_total), eGrid % DEPtotal(eGrid % np_total), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 3')
          DO IP = 1, eGrid % NP_TOTAL
            READ(eGRD % FHNDL, *, IOSTAT = ISTAT) eGrid % XPtotal(IP), eGrid % DEPtotal(IP)
            IF ( ISTAT /= 0 ) CALL WWM_ABORT('error in the grid configuration file')
          END DO
          IF (LVAR1D) THEN
            eGrid % DX1total(0)     = eGrid % XPtotal(2)- eGrid % XPtotal(1)
            eGrid % DX1total(1)     = eGrid % DX1total(0)
            eGrid % DX1total(eGrid%np_total)   = eGrid % XPtotal(eGrid%np_total) - eGrid % XPtotal(eGrid%np_total-1)
            eGrid % DX1total(eGrid%np_total+1) = eGrid % DX1total(eGrid%np_total)
            eGrid % DX2total(0)     = eGrid % DX1total(0)
            eGrid % DX2total(eGrid%np_total+1) = eGrid % DX1total(eGrid%np_total)
            DO IP = 2, eGrid % NP_TOTAL-1 ! Bandwith at gridpoints
              eGrid % DX1total(IP) = (eGrid % XPtotal(IP) - eGrid % XPtotal(IP-1))/2. + (eGrid % XPtotal(IP+1) - eGrid % XPtotal(IP))/2.
            END DO
            DO IP = 2, eGrid % NP_TOTAL ! Stepwidth between gridpoints K and K-1
              eGrid % DX2total(IP) = eGrid % XPtotal(IP) - eGrid % XPtotal(IP-1)
            END DO
            eGrid % DX2total(1) = eGrid % DX1total(0)
          END IF
          CLOSE(eGRD % FHNDL)
        CASE (2)
          IF (IGRIDTYPE == 1) THEN ! system.dat format ... XFN
            OPEN(eGRD % FHNDL, FILE = eGRD % FNAME, STATUS = 'OLD')
            DO I = 1, 2
              READ(eGRD % FHNDL, '(A)') RHEADER
            END DO
            READ(eGRD % FHNDL, *, IOSTAT = ISTAT) ITMP
            IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=1 error in read mnp/mne')
            READ(eGRD % FHNDL, '(A)') RHEADER
            READ(eGRD % FHNDL, *, IOSTAT = ISTAT) JTMP 
            IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=1 error in read mnp/mne')
            EGRID%NP_TOTAL = ITMP + JTMP
            allocate(eGrid % XPtotal(eGrid%np_total), eGrid % YPtotal(eGrid%np_total), eGrid % DEPtotal(eGrid%np_total), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 4')
            DO I = 1, 7
              READ(eGRD % FHNDL, '(A)') RHEADER
            END DO
            DO IP=1,EGRID%NP_TOTAL
              READ(eGRD % FHNDL, *, IOSTAT = ISTAT) KTMP, eGrid % XPtotal(IP), eGrid % YPtotal(IP), eGrid % DEPtotal(IP)
              IF (KTMP+1.ne.IP) THEN
                CALL WWM_ABORT('IGRIDTYPE=1 error in grid read 3')
              ENDIF
              IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=1 error in grid read 4')
            END DO
            DO I = 1, 2
              READ(eGRD % FHNDL, '(A)') RHEADER
            END DO
            READ(eGRD % FHNDL, *, IOSTAT = ISTAT) eGrid % NE_TOTAL
            allocate(eGrid % INEtotal(3, eGrid%ne_total), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 5')
            DO I = 1, 3
              READ(eGRD % FHNDL, '(A)') RHEADER
            END DO
            DO IE=1,EGRID%NE_TOTAL
              READ(eGRD % FHNDL, *, IOSTAT = ISTAT) KTMP, LTMP, MTMP, NTMP, OTMP
              IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=1 error in grid read 5')
              eGrid % INEtotal(1,IE)=KTMP+1
              eGrid % INEtotal(2,IE)=LTMP+1
              eGrid % INEtotal(3,IE)=MTMP+1
            END DO
            CLOSE(eGRD % FHNDL)
          ELSE IF (IGRIDTYPE == 2) THEN ! periodic grid written by mathieu dutour
            OPEN(eGRD%FHNDL, FILE = eGRD%FNAME, STATUS = 'OLD')
            READ(eGRD%FHNDL,*) EGRID%NE_TOTAL, EGRID%NP_TOTAL
            allocate(eGrid % DEPtotal(EGRID%NP_TOTAL), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 6')
            DO IP = 1, EGRID%NP_TOTAL
              READ(eGRD%FHNDL, *, IOSTAT = ISTAT) eGrid % DEPtotal(IP)
              IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=2 error in grid read 1')
            END DO
            allocate(eGrid%TRIAtotal(EGRID%NE_TOTAL), eGrid % INEtotal(3,eGrid%ne_total), eGrid % IENtotal(6,eGrid%ne_total), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 7')
            DO IE = 1, EGRID%NE_TOTAL
              READ(eGRD%FHNDL, *, IOSTAT = ISTAT) eGrid % TRIAtotal(IE)
              IF ( ISTAT /= 0 )  CALL WWM_ABORT('IGRIDTYPE=2 error in grid read 2')
              READ(eGRD%FHNDL, *, IOSTAT = ISTAT) eGrid % INEtotal(:,IE)
              IF ( ISTAT /= 0 )  CALL WWM_ABORT('IGRIDTYPE=2 error in grid read 3')
              READ(eGRD%FHNDL, *, IOSTAT = ISTAT) DXP1, DXP2, DXP3
              IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=2 error in grid read 4')
              READ(eGRD%FHNDL, *, IOSTAT = ISTAT) DYP1, DYP2, DYP3
              IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=2 error in grid read 5')
              eGrid % IENtotal(1,IE) = -DYP2
              eGrid % IENtotal(2,IE) = DXP2
              eGrid % IENtotal(3,IE) = -DYP3
              eGrid % IENtotal(4,IE) = DXP3
              eGrid % IENtotal(5,IE) = -DYP1
              eGrid % IENtotal(6,IE) = DXP1
            END DO
            CLOSE(eGRD%FHNDL)
          ELSE IF (IGRIDTYPE == 3) THEN ! SCHISM gr3
            OPEN(eGRD%FHNDL, FILE = eGRD%FNAME, STATUS = 'OLD')
            READ(eGRD%FHNDL,*)
            READ(eGRD%FHNDL,*, IOSTAT = ISTAT) EGRID%NE_TOTAL, EGRID%NP_TOTAL
            IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=3 error in read mnp/mne')
            allocate(eGrid % XPtotal(eGrid%np_total), eGrid % YPtotal(eGrid%np_total), eGrid % DEPtotal(eGrid%np_total), eGrid % INEtotal(3, eGrid%ne_total), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 8')
            DO IP=1,EGRID%NP_TOTAL
              READ(eGRD%FHNDL, *, IOSTAT = ISTAT) KTMP, XPDTMP, YPDTMP, ZPDTMP
              eGrid % XPtotal(IP)  = XPDTMP
              eGrid % YPtotal(IP)  = YPDTMP
              eGrid % DEPtotal(IP) = ZPDTMP
              IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=3 error in grid reading 1')
            END DO
            DO IE = 1, EGRID%NE_TOTAL
              READ(eGRD%FHNDL, *, IOSTAT = ISTAT) KTMP, LTMP, eGrid % INEtotal(:,IE)
              IF ( ISTAT /= 0 )  CALL WWM_ABORT('IGRIDTYPE=3 error in grid reading 2')
            END DO
            CLOSE(eGRD%FHNDL)
          ELSE IF (IGRIDTYPE == 4) THEN ! Old WWM format
            OPEN(eGRD%FHNDL, FILE = eGRD%FNAME, STATUS = 'OLD')
            READ(eGRD%FHNDL, *, IOSTAT = ISTAT) EGRID%NE_TOTAL 
            READ(eGRD%FHNDL, *, IOSTAT = ISTAT) EGRID%NP_TOTAL 
            allocate(eGrid % XPtotal(eGrid%np_total), eGrid % YPtotal(eGrid%np_total), eGrid % DEPtotal(eGrid%np_total), eGrid%INEtotal(3, eGrid%ne_total), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('allocate error 9')
            DO IP=1,EGRID%NP_TOTAL
              READ(eGRD%FHNDL, *, IOSTAT = ISTAT) eGrid % XPtotal(IP), eGrid % YPtotal(IP), eGrid % DEPtotal(IP)
              IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=4 error in grid read 1')
            END DO
            DO IE=1,EGRID%NE_TOTAL
              READ(eGRD%FHNDL, *, IOSTAT = ISTAT) eGrid % INEtotal(:,IE)
              IF ( ISTAT /= 0 )  CALL WWM_ABORT('IGRIDTYPE=4 error in grid read 2')
            END DO
            CLOSE(eGRD%FHNDL)
#ifdef NCDF
          ELSE IF (IGRIDTYPE == 5) THEN ! Netcdf format
            ISTAT = NF90_OPEN(eGRD%FNAME, NF90_NOWRITE, ncid)
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, ISTAT)

            ISTAT = nf90_inq_varid(ncid, 'ele', var_id)
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, ISTAT)

            ISTAT = nf90_inquire_variable(ncid, var_id, dimids=dimidsB)
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, ISTAT)

            ISTAT = nf90_inquire_dimension(ncid, dimidsB(2), name=MNEstr, len=eGrid%ne_total)
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, ISTAT)
            WRITE(DBG%FHNDL,*) 'EGRID%NE_TOTAL=', EGRID%NE_TOTAL
            WRITE(DBG%FHNDL,*) 'MNEstr=', TRIM(MNEstr)
            FLUSH(DBG%FHNDL)

            ISTAT = nf90_inq_varid(ncid, 'depth', var_id)
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, ISTAT)

            ISTAT = nf90_inquire_variable(ncid, var_id, dimids=dimidsA)
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, ISTAT)

            ISTAT = nf90_inquire_dimension(ncid, dimidsA(1), name=MNPstr, len=eGrid%np_total)
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, ISTAT)
            WRITE(DBG%FHNDL,*) 'EGRID%NP_TOTAL=', EGRID%NP_TOTAL
            WRITE(DBG%FHNDL,*) 'MNPstr=', TRIM(MNPstr)
            FLUSH(DBG%FHNDL)

            allocate(eGrid % XPtotal(eGrid%np_total), eGrid % YPtotal(eGrid%np_total), eGrid % DEPtotal(eGrid%np_total), eGrid%INEtotal(3, eGrid%ne_total), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('allocate error 9')

            ISTAT = nf90_inq_varid(ncid, 'depth', var_id)
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, ISTAT)

            ISTAT = nf90_get_var(ncid, var_id, eGrid % DEPtotal)
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, ISTAT)

            ISTAT = nf90_inq_varid(ncid, 'ele', var_id)
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, ISTAT)

            ISTAT = nf90_get_var(ncid, var_id, eGrid % INEtotal)
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, ISTAT)

            IF (LSPHE) THEN
              ISTAT = nf90_inq_varid(ncid, 'lon', var_id1)
              CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, ISTAT)

              ISTAT = nf90_inq_varid(ncid, 'lat', var_id2)
              CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, ISTAT)
            ELSE
              ISTAT = nf90_inq_varid(ncid, 'x', var_id1)
              CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, ISTAT)

              ISTAT = nf90_inq_varid(ncid, 'y', var_id2)
              CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 15, ISTAT)
            END IF
            ISTAT = nf90_get_var(ncid, var_id1, eGrid % XPtotal)
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 16, ISTAT)

            ISTAT = nf90_get_var(ncid, var_id2, eGrid % YPtotal)
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 17, ISTAT)

            ISTAT = NF90_CLOSE(ncid)
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 18, ISTAT)
#endif
          ELSE
            CALL WWM_ABORT('IGRIDTYPE WRONG')
          END IF
          DO IE=1,EGRID%NE_TOTAL
            DO I=1,3
              IP=eGrid % INEtotal(I,IE)
              IF ((IP .lt. 1).or.(IP.gt.EGRID%NP_TOTAL)) THEN
                Print *, 'IE=', IE, ' I=', I
                Print *, 'IP=', IP, ' EGRID%NP_TOTAL=', EGRID%NP_TOTAL
                CALL WWM_ABORT('INCOHERENCY IN INEtotal')
              END IF
            END DO
          END DO
        CASE DEFAULT
          CALL WWM_ABORT('WRONG GRID DIMENSION')
      END SELECT
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_SPATIAL_GRID_TOTAL_KERNEL(eGrid, DimMode, LVAR1D, Lsphe, eGRD, iGridType)
      USE DATAPOOL, only : rkind, istat, GridInformation, FILEDEF
      IMPLICIT NONE
      type(GridInformation), intent(out) :: eGrid
      integer, intent(in) :: DimMode
      logical, intent(in) :: LVAR1D
      logical, intent(in) :: Lsphe
      type(FILEDEF), intent(in) :: eGRD
      integer, intent(in) :: iGridType
      !
      integer :: rbuf_int(2)
      real(rkind), allocatable :: rbuf_real(:)
      integer iProc, IP, IE, nb_real, idx
#ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_GRID) THEN
        CALL SINGLE_READ_SPATIAL_GRID_TOTAL(eGrid, DimMode, LVAR1D, Lsphe, eGRD, iGridType)
      ELSE
        IF (DIMMODE .ne. 2) THEN
          CALL WWM_ABORT('Parallel mode only for 2D')
        ENDIF
        IF (myrank .eq. 0) THEN
          CALL SINGLE_READ_SPATIAL_GRID_TOTAL(eGrid, DimMode, LVAR1D, Lsphe, eGRD, iGridType)
          rbuf_int(1)=eGrid % np_total
          rbuf_int(2)=eGrid % ne_total
          DO iProc=2,nproc
            CALL MPI_SEND(rbuf_int,2,itype, iProc-1, 30, comm, ierr)
          END DO
          DO iProc=2,nproc
            CALL MPI_SEND(INEtotal,3*eGrid % ne_total,itype, iProc-1, 32, comm, ierr)
          END DO
          IF (IGRIDTYPE .eq. 2) THEN
            nb_real=eGrid % np_total + 7*eGrid % ne_total
            allocate(rbuf_real(nb_real), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('allocate error 10')
            idx=0
            DO IP=1,NP_TOTAL
              idx=idx+1
              rbuf_real(idx)=eGrid % DEPtotal(IP)
            END DO
            DO IE=1,eGrid % NE_TOTAL
              rbuf_real(idx+1)=eGrid % TRIAtotal(IE)
              rbuf_real(idx+2)=eGrid % IENtotal(1,IE)
              rbuf_real(idx+3)=eGrid % IENtotal(2,IE)
              rbuf_real(idx+4)=eGrid % IENtotal(3,IE)
              rbuf_real(idx+5)=eGrid % IENtotal(4,IE)
              rbuf_real(idx+6)=eGrid % IENtotal(5,IE)
              rbuf_real(idx+7)=eGrid % IENtotal(6,IE)
              idx=idx+7
            END DO
          ELSE
            nb_real=3*eGrid % np_total
            allocate(rbuf_real(nb_real), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('allocate error 11')
            idx=0
            DO IP=1,NP_TOTAL
              rbuf_real(idx+1)=eGrid % XPtotal(IP)
              rbuf_real(idx+2)=eGrid % YPtotal(IP)
              rbuf_real(idx+3)=eGrid % DEPtotal(IP)
              idx=idx+3
            END DO
          END IF
          DO iProc=2,nproc
            CALL MPI_SEND(rbuf_real,nb_real,rtype, iProc-1, 34, comm, ierr)
          END DO
          deallocate(rbuf_real)
        ELSE
          CALL MPI_RECV(rbuf_int,2,itype, 0, 30, comm, istatus, ierr)
          eGrid % np_total=rbuf_int(1)
          eGrid % ne_total=rbuf_int(2)
          allocate(eGrid % INEtotal(3,eGrid % ne_total), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('allocate error 12')
          CALL MPI_RECV(eGrid % INEtotal,3*eGrid % ne_total,itype, 0, 32, comm, istatus, ierr)
          IF (IGRIDTYPE .eq. 2) THEN
            nb_real=eGrid % np_total + 7*eGrid % ne_total
            allocate(rbuf_real(nb_real), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('allocate error 13')
            CALL MPI_RECV(rbuf_real,nb_real,rtype, 0, 34, comm, istatus, ierr)
            allocate(eGrid % DEPtotal(eGrid % np_total), eGrid % TRIAtotal(eGrid % ne_total), eGrid % IENtotal(6,eGrid % ne_total), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('allocate error 14')
            idx=0
            DO IP=1,eGrid % NP_TOTAL
              idx=idx+1
              eGrid % DEPtotal(IP)=rbuf_real(idx)
            END DO
            DO IE=1,eGrid % NE_TOTAL
              eGrid % TRIAtotal(IE)=rbuf_real(idx+1)
              eGrid % IENtotal(1,IE)=rbuf_real(idx+2)
              eGrid % IENtotal(2,IE)=rbuf_real(idx+3)
              eGrid % IENtotal(3,IE)=rbuf_real(idx+4)
              eGrid % IENtotal(4,IE)=rbuf_real(idx+5)
              eGrid % IENtotal(5,IE)=rbuf_real(idx+6)
              eGrid % IENtotal(6,IE)=rbuf_real(idx+7)
              idx=idx+7
            END DO
          ELSE
            nb_real=3*np_total
            allocate(rbuf_real(nb_real), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('allocate error 15')
            allocate(eGrid % DEPtotal(np_total), eGrid % XPtotal(np_total), eGrid % YPtotal(np_total), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('allocate error 16')
            CALL MPI_RECV(rbuf_real,nb_real,rtype, 0, 34, comm, istatus, ierr)
            idx=0
            DO IP=1,NP_TOTAL
              eGrid % XPtotal(IP) = rbuf_real(idx+1)
              eGrid % YPtotal(IP) = rbuf_real(idx+2)
              eGrid % DEPtotal(IP)= rbuf_real(idx+3)
              idx=idx+3
            END DO
          END IF
          deallocate(rbuf_real)
        END IF
      END IF
#else
      CALL SINGLE_READ_SPATIAL_GRID_TOTAL(eGrid, DimMode, LVAR1D, Lsphe, eGRD, iGridType)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_SPATIAL_GRID_TOTAL
      USE DATAPOOL
      IMPLICIT NONE
      type(GridInformation) eGrid
      integer IP, IE
      CALL READ_SPATIAL_GRID_TOTAL_KERNEL(eGrid, DIMMODE, LVAR1D, LSPHE, GRD, IGRIDTYPE)
      np_total = eGrid % np_total
      ne_total = eGrid % ne_total
      IF (DIMMODE .eq. 1) THEN
        allocate(XPtotal(np_total), DEPtotal(np_total), DX1total(0:np_total+1), DX2total(0:np_total+1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('allocate error')
        DO IP=1,np_total
          XPtotal(IP)  = eGrid % XPtotal(IP)
          DEPtotal(IP) = eGrid % DEPtotal(IP)
        END DO
        DO IP=0,np_total+1
          DX1total(IP) = eGrid % DX1total(IP)
          DX2total(IP) = eGrid % DX2total(IP)
        END DO
      ELSE
        IF (IGRIDTYPE .eq. 2) THEN
          allocate(DEPtotal(NP_TOTAL), TRIAtotal(NE_TOTAL), INEtotal(3,ne_total), IENtotal(6,ne_total), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 6')
          DO IP=1,np_total
            DEPtotal(IP) = eGrid % DEPtotal(IP)
          END DO
          DO IE=1,ne_total
            TRIAtotal(IE) = eGrid % TRIAtotal(IE)
            INEtotal(:,IE) = eGrid % INEtotal(:,IE)
            IENtotal(:,IE) = eGrid % IENtotal(:,IE)
          END DO
        ELSE
          allocate(XPtotal(np_total), YPtotal(np_total), DEPtotal(np_total), INEtotal(3, eGrid%ne_total), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 4')
          DO IP=1,np_total
            XPtotal(IP) = eGrid % XPtotal(IP)
            YPtotal(IP) = eGrid % YPtotal(IP)
            DEPtotal(IP) = eGrid % DEPtotal(IP)
          END DO
          DO IE=1,ne_total
            INEtotal(:,IE) = eGrid % INEtotal(:,IE)
          END DO
        END IF
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXPORT_GRID_WW3_FORMAT
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER nbDirichlet, nbIsland
      INTEGER :: MapDirect(np_total)
      INTEGER :: MapRevert(np_total)
      INTEGER, allocatable :: IPbound(:), IPisland(:), ACTIVE(:)
      INTEGER IP, IE2, eVal, IE, idxDirichlet, idxIsland
#ifdef MPI_PARALL_GRID      
      IF (myrank .ne. 0) THEN
        RETURN
      END IF
#endif
      nbDirichlet=0
      nbIsland=0
      DO IP=1,np_total
        eVal=IOBPtotal(IP)
        IF (eVal .eq. 2) THEN
          nbDirichlet=nbDirichlet+1
        END IF
        IF ((eVal .eq. 1).or.(eVal .eq. 3).or.(eVal .eq. 4)) THEN
          nbIsland=nbIsland+1
        END IF
      END DO
      allocate(IPbound(nbDirichlet), IPisland(nbIsland), ACTIVE(nbDirichlet), stat=istat)
      ACTIVE(:)=1    ! right now we are proceding this way. Maybe there is work here to improve
      IF (istat/=0) CALL WWM_ABORT('allocate error 17')
      idxDirichlet=0
      idxIsland=0
      DO IP=1,np_total
        eVal=IOBPtotal(IP)
        IF (eVal .eq. 2) THEN
          idxDirichlet=idxDirichlet+1
          IPbound(idxDirichlet)=IP
        END IF
        IF ((eVal .eq. 1).or.(eVal .eq. 3).or.(eVal .eq. 4)) THEN
          idxIsland=idxIsland+1
          IPisland(idxIsland)=IP
        END IF
      END DO
      OPEN(FHNDL_EXPORT_GRID_WW3, FILE = 'mesh.msh', STATUS='unknown')
      !
      !*** write gmsh header
      !
      WRITE(FHNDL_EXPORT_GRID_WW3,'(A)') '$MeshFormat'
      WRITE(FHNDL_EXPORT_GRID_WW3,'(A)') '2 0 8'
      WRITE(FHNDL_EXPORT_GRID_WW3,'(A)') '$EndMeshFormat'
      WRITE(FHNDL_EXPORT_GRID_WW3,'(A)') '$Nodes'
      WRITE(FHNDL_EXPORT_GRID_WW3,*) np_total
      DO IP = 1, np_total
         WRITE(FHNDL_EXPORT_GRID_WW3,'(I10,3F30.20)') IP, XPtotal(IP), YPtotal(IP), DEPtotal(IP)
      END DO
      WRITE(FHNDL_EXPORT_GRID_WW3,'(A)') '$EndNodes'
      WRITE(FHNDL_EXPORT_GRID_WW3,'(A)') '$Elements'
      WRITE(FHNDL_EXPORT_GRID_WW3,*)  ne_total + nbDirichlet + nbIsland
      !
      !*** write boundary
      !
      IE2=1
      DO IP = 1, nbDirichlet
        WRITE(FHNDL_EXPORT_GRID_WW3,'(6I10)') IE2, 15, 2, ACTIVE(IP), 0, IPbound(IP)
        IE2=IE2+1
      END DO
      !
      !*** write island
      !
      DO IP = 1, nbIsland
        WRITE(FHNDL_EXPORT_GRID_WW3,'(6I10)') IE2, 15, 2, 0, IP, IPisland(IP)
        IE2=IE2+1
      END DO
      !
      !*** write gmsh elements
      !
      DO IE = 1, ne_total
         WRITE(FHNDL_EXPORT_GRID_WW3,'(9I8)') IE2, 2, 3, 0, IE, 0, INE(1,IE), IEN(2,IE), IEN(3,IE)
         IE2=IE2+1
      END DO
      WRITE(FHNDL_EXPORT_GRID_WW3,'(A)') '$EndElements'
      CLOSE(FHNDL_EXPORT_GRID_WW3)
      deallocate(IPbound, IPisland, ACTIVE)
      END SUBROUTINE
      
