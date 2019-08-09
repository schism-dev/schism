#include "wwm_functions.h"
#ifdef NCDF
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_NESTING
      USE DATAPOOL
      IMPLICIT NONE
      real(rkind) :: XYTMP(2,MNP)
      real(rkind) :: eWI(3)
      integer iGrid, np_total_loc, IWBMNP_loc
      type(FILEDEF) eGRD, eBND
      integer IP, idx
      REAL(rkind) eX, eY
      integer eElt, NI(3), eIdx
      !
      ! First reading the grids
      !
      ALLOCATE(ListNestInfo(NB_GRID_NEST), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 1')
      DO iGrid=1,NB_GRID_NEST
        IGRIDTYPE = ListIGRIDTYPE(iGrid)
        eGRD % FNAME = ListFILEGRID(iGrid)
        eGRD % FHNDL = 24037
        CALL READ_SPATIAL_GRID_TOTAL_KERNEL(ListNestInfo(iGrid) % eGrid, DIMMODE, LVAR1D, LSPHE, eGRD, IGRIDTYPE)
        !
!        np_total_loc = ListNestInfo(iGrid) % eGrid % np_total
!        allocate(ListNestInfo(iGrid) % NodeBelonging(np_total_loc), stat=istat)
!        IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 2')
        !
        allocate(ListNestInfo(iGrid) % IOBPtotal(np_total_loc), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 3')
        eBND % FNAME = ListFILEBOUND(iGrid)
        eBND % FHNDL = 24977
        CALL SINGLE_READ_IOBP_TOTAL(IOBPtotal, IGRIDTYPE, eBND, np_total_loc)
        !
        IWBMNP_loc=0
        DO IP=1,np_total_loc
          IF ((IOBPtotal(IP) .eq. 2) .or. (IOBPtotal(IP) .eq. 3)) THEN
            IWBMNP_loc = IWBMNP_loc + 1
          END IF
        END DO
        ListNestInfo(iGrid) % IWBMNP = IWBMNP_loc
        !
        allocate(ListNestInfo(iGrid) % IWBNDLC(IWBMNP_loc), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 2')
        idx=0
        DO IP=1,np_total_loc
          IF ((IOBPtotal(IP) .eq. 2) .or. (IOBPtotal(IP) .eq. 3)) THEN
            idx=idx+1
            ListNestInfo(iGrid) % IWBNDLC(idx) = IP
          END IF
        END DO
        !
        
      END DO
      !
      ! Now we find the IE and weights for interpolation
      !
      XYTMP(1,:) = XP
      XYTMP(2,:) = YP
      DO iGrid=1,NB_GRID_NEST
        np_total_loc = ListNestInfo(iGrid) % eGrid % np_total
        IWBMNP_loc = ListNestInfo(iGrid) % IWBMNP
        IF (L_HOTFILE) THEN
          allocate(ListNestInfo(iGrid) % HOT_IE(np_total_loc),ListNestInfo(iGrid) % HOT_W(3, np_total_loc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 2')
          DO IP=1,np_total_loc
            eX=ListNestInfo(iGrid) % eGrid % XPtotal(IP)
            eY=ListNestInfo(iGrid) % eGrid % YPtotal(IP)
            CALL FIND_ELE(MNE,MNP,INE,XYTMP,eX, eY, eElt)
            ListNestInfo(iGrid) % HOT_IE(IP) = eElt
            IF (eElt .gt. 0) THEN
              NI                 = INE(:,eElt)
              CALL INTELEMENT_COEF(XP(NI),YP(NI), eX, eY, eWI)
              ListNestInfo(iGrid) % HOT_W(:, IP) = eWI
            END IF
          END DO
        END IF
        IF (L_BOUC_PARAM .or. L_BOUC_SPEC) THEN
          allocate(ListNestInfo(iGrid) % BOUC_IE(IWBMNP_loc),ListNestInfo(iGrid) % BOUC_W(3, IWBMNP_loc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 2')
          DO eIdx=1,IWBMNP_loc
            IP=ListNestInfo(iGrid) % IWBNDLC(eIdx)
            eX=ListNestInfo(iGrid) % eGrid % XPtotal(IP)
            eY=ListNestInfo(iGrid) % eGrid % YPtotal(IP)
            CALL FIND_ELE(MNE, MNP, INE, XYTMP, eX, eY, eElt)
            ListNestInfo(iGrid) % BOUC_IE(IP) = eElt
            IF (eElt .gt. 0) THEN
              NI                 = INE(:,eElt)
              CALL INTELEMENT_COEF(XP(NI),YP(NI), eX, eY, eWI)
              ListNestInfo(iGrid) % BOUC_W(:, IP) = eWI
            END IF
          END DO
        END IF
      END DO
      !
      ! Now we timings needed by the model
      !
      DO iGRid=1,NB_GRID_NEST
        ListNestInfo(iGrid) % eTime % BEGT = ListBEGTC(iGrid)
        ListNestInfo(iGrid) % eTime % DELT = ListDELTC(iGrid)
        ListNestInfo(iGrid) % eTime % UNIT = ListUNITC(iGrid)
        ListNestInfo(iGrid) % eTime % ENDT = ListENDTC(iGrid)
        CALL CT2MJD(ListNestInfo(iGrid) % eTime % BEGT, ListNestInfo(iGrid) % eTime % BMJD)
        CALL CT2MJD(ListNestInfo(iGrid) % eTime % ENDT, ListNestInfo(iGrid) % eTime % EMJD)
        CALL CU2SEC(ListNestInfo(iGrid) % eTime % UNIT, ListNestInfo(iGrid) % eTime % DELT)

        ListNestInfo(iGrid) % eTime % TOTL = (ListNestInfo(iGrid) % eTime % EMJD - ListNestInfo(iGrid) % eTime % BMJD) * DAY2SEC
        ListNestInfo(iGrid) % eTime % ISTP = NINT(ListNestInfo(iGrid) % eTime % TOTL / ListNestInfo(iGrid) % eTime % DELT) + 1
        ListNestInfo(iGrid) % eTime % TMJD = ListNestInfo(iGrid) % eTime % BMJD
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NESTING_OUTPUT_HOTFILE(iGrid)
      USE DATAPOOL
      USE WWM_HOTFILE_MOD
      IMPLICIT NONE
      integer, intent(in) :: iGrid
      character(len=140) FILERET
      real(rkind), allocatable :: ACwrite(:,:,:), VAR_ONEDwrite(:,:)
      real(rkind), allocatable :: ACsend(:,:,:), VAR_ONEDsend(:,:)
      integer, allocatable :: ListStatus(:)
      real(rkind) :: eVect(nbOned)
      !
      integer eInt(1)
      integer np_write, ne_write
      integer MULTIPLEOUT_W
      logical GRIDWRITE_W, IOBPD_HISTORY_W, WriteOutputProcess
      integer nbMatch, nbTime
      integer IP, IE, I, idx, IP2, iProc, nbMatchLoc
      real(rkind) eW
      integer, allocatable :: ListMatch(:)
      real(rkind) eTimeDay
      integer POS
      np_write=ListNestInfo(iGrid) % eGrid % np_total
      ne_write=ListNestInfo(iGrid) % eGrid % ne_total
#ifdef MPI_PARALL_GRID
      IF (myrank .eq. 0) THEN
#endif
        FILERET = ListPrefix(iGrid) // '_hotfile.nc'
        MULTIPLEOUT_W = .FALSE.
        GRIDWRITE_W = .FALSE.
        IOBPD_HISTORY_W = .FALSE.
        WriteOutputProcess = .TRUE.
        nbTime=-1
        CALL WRITE_HOTFILE_PART_1(FILERET, nbTime, MULTIPLEOUT_W, GRIDWRITE_W, IOBPD_HISTORY_W, np_write, ne_write)
        CALL WRITE_NETCDF_HEADERS_2(FILERET, MULTIPLEOUT_W, WriteOutputProcess, GRIDWRITE_W, np_write, ne_write)
        allocate(ACwrite(MSC,MDC,np_write), VAR_ONEDwrite(nbOned, np_write), ListStatus(np_write), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 1')
        ListStatus=0
#ifdef MPI_PARALL_GRID
      END IF
#endif
      nbMatch=0
      DO IP=1,np_write
        IE=ListNestInfo(iGrid) % HOT_IE(IP)
        IF (IE .gt. 0) THEN
          nbMatch = nbMatch + 1
        END IF
      END DO
      allocate(Listmatch(nbMatch), ACsend(MSC,MDC,nbMatch), VAR_ONEDsend(nbOned, nbMatch), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 1')
      ACsend=0
      VAR_ONEDsend=0
      idx=0
      DO IP=1,np_write
        IE=ListNestInfo(iGrid) % HOT_IE(IP)
        IF (IE .gt. 0) THEN
          idx=idx+1
          ListMatch(idx)=IP
          DO I=1,3
            eW=ListNestInfo(iGrid) % HOT_W(I,IP)
            IP2=INE(I,IE)
            eVect(1)=WINDXY(IP2,1)
            eVect(2)=WINDXY(IP2,2)
            eVect(3)=PRESSURE(IP2)
            eVect(4)=DVWIND(IP2,1)
            eVect(5)=DVWIND(IP2,2)
            eVect(6)=CURTXY(IP2,1)
            eVect(7)=CURTXY(IP2,2)
            eVect(8)=DVCURT(IP2,1)
            eVect(9)=DVCURT(IP2,2)
            eVect(10)=DDEP(IP2,1)
            eVect(11)=DDEP(IP2,2)
            eVect(12)=DCUX(IP2,1)
            eVect(13)=DCUX(IP2,2)
            eVect(14)=DCUY(IP2,1)
            eVect(15)=DCUY(IP2,2)
            eVect(16)=WATLEV(IP2)
            eVect(17)=WATLEVOLD(IP2)
            eVect(18)=WLDEP(IP2)
            eVect(19)=DEPDT(IP2)
            eVect(20)=QBLOCAL(IP2)
            eVect(21)=DISSIPATION(IP2)
            eVect(22)=AIRMOMENTUM(IP2)
            eVect(23)=UFRIC(IP2)
            eVect(24)=ALPHA_CH(IP2)
            eVect(25)=TAUW(IP2)
            eVect(26)=TAUTOT(IP2)
            eVect(27)=TAUWX(IP2)
            eVect(28)=TAUWY(IP2)
            eVect(29)=TAUHF(IP2)
            eVect(30)=Z0(IP2)
            eVect(31)=CD(IP2)
            eVect(32)=USTDIR(IP2)
            eVect(33)=RSXX(IP2)
            eVect(34)=RSXY(IP2)
            eVect(35)=RSYY(IP2)
            eVect(36)=FORCEXY(IP2,1)
            eVect(37)=FORCEXY(IP2,2)
            ACsend(:,:,idx) = ACsend(:,:,idx) + eW * AC2(:,:,IP2)
            VAR_ONEDsend(:,idx) = VAR_ONEDsend(:,idx) + eW * eVect
          END DO
        END IF
      END DO
#ifdef MPI_PARALL_GRID
      IF (myrank .eq. 0) THEN
        DO idx=1,nbMatch
          IP=ListMatch(idx)
          ACwrite(:,:,IP)     = ACsend(:,:,idx)
          VAR_ONEDwrite(:,IP) = VAR_ONEDsend(:,idx)
          ListStatus(IP)=1
        END DO
        deallocate(ACsend, VAR_ONEDsend, ListMatch)
        DO iProc=2,nproc
          CALL MPI_RECV(eInt, 1, itype, iProc-1, 2401, comm, istatus, ierr)
          nbMatchLoc=eInt(1)
          allocate(ListMatch(nbMatchLoc), ACsend(MSC,MDC,nbMatchLoc), VAR_ONEDsend(nbOned, nbMatchLoc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 1')
          CALL MPI_RECV(ListMatch, nbMatchLoc, itype, iProc-1, 2402, comm, istatus, ierr)
          CALL MPI_RECV(ACsend, MSC*MDC*nbMatchLoc, rtype, iProc-1, 2403, comm, istatus, ierr)
          CALL MPI_RECV(VAR_ONEDsend, nbOned*nbMatchLoc, rtype, iProc-1, 2404, comm, istatus, ierr)
          DO idx=1,nbMatchLoc
            IP=ListMatch(idx)
            ACwrite(:,:,IP) = ACsend(:,:,idx)
            VAR_ONEDwrite(:,IP) = VAR_ONEDsend(:,idx)
            ListStatus(IP)=1
          END DO
          deallocate(ListMatch, ACsend, VAR_ONEDsend)
        END DO
      ELSE
        eInt(1)=nbMatch
        CALL MPI_SEND(eInt, 1, itype, 0, 2401, comm, ierr)
        CALL MPI_SEND(ListMatch, nbMatch, itype, 0, 2402, comm, ierr)
        CALL MPI_SEND(ACsend, MSC*MDC*nbMatch, rtype, 0, 2403, comm, ierr)
        CALL MPI_SEND(VAR_ONEDsend, nbOned*nbMatch, rtype, 0, 2404, comm, ierr)
      END IF
#else
      ACwrite = ACsend
      VAR_ONEDwrite = VAR_ONEDsend
#endif      
#ifdef MPI_PARALL_GRID
      IF (myrank .eq. 0) THEN
#endif
        eTimeDay=ListNestInfo(iGrid) % eTime % BMJD
        POS=1
        CALL WRITE_HOTFILE_PART_2(FILERET, eTimeDay, POS, np_write, ACwrite, VAR_ONEDwrite)
#ifdef MPI_PARALL_GRID
      END IF
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NESTING_BOUNDARY_CONDITION(iGrid)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
      integer, intent(in) :: iGrid
      character (len = *), parameter :: CallFct = "NESTING_BOUNDARY_CONDITION"
      real(rkind), allocatable :: WBACwrite(:,:,:), SPPARMwrite(:,:)
      real(rkind), allocatable :: WBACsend(:,:,:), SPPARMsend(:,:)
      integer, allocatable :: ListStatus(:)
      REAL(rkind) :: CURTXYLOC(2), DEPLOC, WATLEVLOC, WKLOC(MSC), ACLOC(MSC,MDC)
      real(rkind) :: eVect(8)
      real(rkind) :: WVK,WVCG,WVKDEP,WVN,WVC
      integer IE, nbTime, IP, idx
      integer nbMatch, np_write, nbBound
      integer, allocatable :: ListMatch(:)
      integer IP2, I, IS
      integer ISMAX, nbMatchLoc
      integer recs_his, irec_dim, var_id
      real(rkind) HS, TM01, TM10, TM02, KLM, WLM
      real(rkind) ETOTS, ETOTC, DM, DSPR
      real(rkind) eW
      integer eInt(1), iProc
      integer iret, ncid
      real(rkind) eTimeDay
      character(len=140) FILERET
      np_write=ListNestInfo(iGrid) % eGrid % np_total
      nbBound=ListNestInfo(iGrid) % IWBMNP
#ifdef MPI_PARALL_GRID
      IF (myrank .eq. 0) THEN
#endif
        FILERET = ListPrefix(iGrid) // '_boundary.nc'
        nbTime=-1 
        CALL WRITE_NETCDF_BOUND_HEADERS_1(FILERET, nbTime, np_write, nbBound, L_BOUC_PARAM, L_BOUC_SPEC)
        CALL WRITE_NETCDF_BOUND_HEADERS_2(FILERET, np_write, ListNestInfo(iGrid) % IOBPtotal, nbBound, ListNestInfo(iGrid) % IWBNDLC)
        IF (L_BOUC_PARAM) THEN
          allocate(SPPARMwrite(8,nbBound), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 1')
        END IF
        IF (L_BOUC_SPEC) THEN
          allocate(WBACwrite(MSC,MDC,nbBound), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 1')
        END IF
        allocate(ListStatus(nbBound), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 1')
        ListStatus=0
#ifdef MPI_PARALL_GRID
      END IF
#endif
      nbMatch=0
      DO IP=1,nbBound
        IE=ListNestInfo(iGrid) % BOUC_IE(IP)
        IF (IE .gt. 0) THEN
          nbMatch = nbMatch + 1
        END IF
      END DO
      IF (L_BOUC_PARAM) THEN
        allocate(SPPARMsend(8,nbMatch), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 1')
        SPPARMsend=0
      END IF
      IF (L_BOUC_SPEC) THEN
        allocate(WBACsend(MSC,MDC,nbMatch), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 1')
        WBACsend=0
      END IF
      allocate(Listmatch(nbMatch), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 1')
      idx=0
      DO IP=1,nbBound
        IE=ListNestInfo(iGrid) % BOUC_IE(IP)
        IF (IE .gt. 0) THEN
          idx=idx+1
          ListMatch(idx)=IP
          ACLOC=0
          DEPLOC=0
          DO I=1,3
            eW=ListNestInfo(iGrid) % HOT_W(I,IP)
            IP2=INE(I,IE)
            ACLOC = ACLOC + eW * AC2(:,:,IP2)
            IF (L_BOUC_PARAM) THEN
              DEPLOC = DEPLOC + eW * DEP(IP2)
              CURTXYLOC = CURTXYLOC + eW * CURTXY(IP2,:)
              WATLEVLOC = WATLEVLOC + eW * WATLEV(IP2)
            END IF
          END DO
          IF (L_BOUC_SPEC) THEN
            WBACsend(:,:,idx) = ACLOC
          END IF
          IF (L_BOUC_PARAM) THEN
            DO IS = 1, MSC
              CALL ALL_FROM_TABLE(SPSIG(IS),DEPLOC,WVK,WVCG,WVKDEP,WVN,WVC)
              WKLOC(IS) = WVK
            END DO
            ISMAX=MSC
            CALL MEAN_PARAMETER_LOC(ACLOC,CURTXYLOC,DEPLOC,WKLOC,ISMAX,HS,TM01,TM02,TM10,KLM,WLM)
            CALL MEAN_DIRECTION_AND_SPREAD_LOC(ACLOC,ISMAX,ETOTS,ETOTC,DM,DSPR)
            eVect(1) = HS
            eVect(2) = TM01
            eVect(3) = DM
            eVect(4) = DSPR
            eVect(5) = 2  ! JONSWAP
            eVect(6) = 2 ! directional spreading in degrees
            eVect(7) = 0.1 ! not used
            eVect(8) = 3.3 ! peak enhancement factor of JONSWAP
            SPPARMsend(:,idx) = eVect
          END IF
        END IF
      END DO
#ifdef MPI_PARALL_GRID
      IF (myrank .eq. 0) THEN
        DO idx=1,nbMatch
          IP=ListMatch(idx)
          IF (L_BOUC_SPEC) THEN
            WBACwrite(:,:,IP) = WBACsend(:,:,idx)
          END IF
          IF (L_BOUC_PARAM) THEN
            SPPARMwrite(:,IP) = SPPARMsend(:,idx)
          END IF
          ListStatus(IP)=1
        END DO
        IF (L_BOUC_SPEC) THEN
          deallocate(WBACsend)
        END IF
        IF (L_BOUC_PARAM) THEN
          deallocate(SPPARMsend)
        END IF
        deallocate(ListMatch)
        DO iProc=2,nproc
          CALL MPI_RECV(eInt, 1, itype, iProc-1, 2401, comm, istatus, ierr)
          nbMatchLoc=eInt(1)
          allocate(ListMatch(nbMatchLoc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 1')
          CALL MPI_RECV(ListMatch, nbMatchLoc, itype, iProc-1, 2402, comm, istatus, ierr)
          IF (L_BOUC_SPEC) THEN
            allocate(WBACsend(MSC,MDC,nbMatchLoc), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 1')
            CALL MPI_RECV(WBACsend, MSC*MDC*nbMatchLoc, rtype, iProc-1, 2403, comm, istatus, ierr)
            DO idx=1,nbMatchLoc
              IP=ListMatch(idx)
              WBACwrite(:,:,IP) = WBACsend(:,:,idx)
              ListStatus(IP)=1
            END DO
            deallocate(WBACsend)
          END IF
          IF (L_BOUC_PARAM) THEN
            allocate(SPPARMsend(8, nbMatchLoc), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 1')
            CALL MPI_RECV(SPPARMsend, 8*nbMatchLoc, rtype, iProc-1, 2404, comm, istatus, ierr)
            DO idx=1,nbMatchLoc
              IP=ListMatch(idx)
              SPPARMwrite(:,IP) = SPPARMsend(:,idx)
              ListStatus(IP)=1
            END DO
            deallocate(SPPARMsend)
          END IF
          deallocate(ListMatch)
        END DO
      ELSE
        eInt(1)=nbMatch
        CALL MPI_SEND(eInt, 1, itype, 0, 2401, comm, ierr)
        CALL MPI_SEND(ListMatch, nbMatch, itype, 0, 2402, comm, ierr)
        IF (L_BOUC_SPEC) THEN
          CALL MPI_SEND(WBACsend, MSC*MDC*nbMatch, rtype, 0, 2403, comm, ierr)
        END IF
        IF (L_BOUC_PARAM) THEN
          CALL MPI_SEND(SPPARMsend, 8*nbMatch, rtype, 0, 2404, comm, ierr)
        END IF
      END IF
#else
      WBACwrite = WBACsend
      SPPARMwrite = SPPARMsend
#endif      
#ifdef MPI_PARALL_GRID
      IF (myrank .eq. 0) THEN
#endif
        eTimeDay = ListNestInfo(iGrid) % eTime % TMJD
        iret=nf90_open(TRIM(FILERET), NF90_WRITE, ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
        iret=nf90_inquire(ncid, unlimitedDimId = irec_dim)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
        iret=nf90_inquire_dimension(ncid, irec_dim,len = recs_his)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, iret)
        recs_his=recs_his+1
        CALL WRITE_NETCDF_TIME(ncid, recs_his, eTimeDay)
        IF (L_BOUC_PARAM) THEN
          iret=nf90_inq_varid(ncid, 'SPPARM', var_id)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, iret)
          IF (NF90_RUNTYPE == NF90_OUTTYPE_BOUC) THEN
            iret=nf90_put_var(ncid,var_id,SPPARMwrite, start=(/1,1,recs_his/), count = (/8, nbBound,1/))
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, iret)
          ELSE
            iret=nf90_put_var(ncid,var_id,SNGL(SPPARMwrite), start=(/1,1,recs_his/), count = (/8, nbBound,1/))
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 10, iret)
          ENDIF
        END IF
        IF (L_BOUC_SPEC) THEN
          !
          iret=nf90_inq_varid(ncid, 'WBAC', var_id)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 11, iret)
          IF (NF90_RUNTYPE == NF90_OUTTYPE_BOUC) THEN
            iret=nf90_put_var(ncid,var_id,WBACwrite, start=(/1,1,1,recs_his/), count = (/MSC,MDC, nbBound,1/))
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 12, iret)
          ELSE
            iret=nf90_put_var(ncid,var_id,SNGL(WBACwrite), start=(/1,1,1,recs_his/), count = (/MSC,MDC, nbBound,1/))
            CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 13, iret)
          ENDIF
        END IF
        iret=nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 14, iret)
        IF (L_BOUC_PARAM) THEN
          deallocate(SPPARMwrite)
        END IF
        IF (L_BOUC_SPEC) THEN
          deallocate(WBACwrite)
        END IF
#ifdef MPI_PARALL_GRID
      END IF
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DO_NESTING_OPERATION
      USE DATAPOOL
      IMPLICIT NONE
      LOGICAL, SAVE :: INITDONE = .FALSE.
      integer iGrid
      real(rkind) DeltaTimeDiff
      !
      ! First the init
      !
      IF (INITDONE .eqv. .FALSE.) THEN
        CALL INIT_NESTING
      END IF
      INITDONE = .TRUE.
      !
      ! Then the HOTFILE
      !
      DO iGrid=1,NB_GRID_NEST
        IF ((MAIN%TMJD .GE. ListNestInfo(iGrid) % eTime % TMJD - 1.E-8) .AND. (MAIN%TMJD .LE. ListNestInfo(iGrid) % eTime % EMJD)) THEN
          DeltaTimeDiff = abs(MAIN % TMJD - ListNestInfo(iGrid) % eTime % BMJD)
          IF (L_HOTFILE .and. DeltaTimeDiff .le. 1.e-8) THEN
            CALL NESTING_OUTPUT_HOTFILE(iGrid)
          END IF
          IF (L_BOUC_PARAM .or. L_BOUC_SPEC) THEN
            CALL NESTING_BOUNDARY_CONDITION(iGrid)
          END IF
          ListNestInfo(iGrid) % eTime % TMJD = ListNestInfo(iGrid) % eTime % TMJD + ListNestInfo(iGrid) % eTime % DELT*SEC2DAY
        END IF
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#endif
