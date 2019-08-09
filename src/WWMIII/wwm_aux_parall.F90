#include "wwm_functions.h"
!**********************************************************************
!* Determine if nodes are in the domain for the computation of        *
!* energy and similar things                                          *
!**********************************************************************
      SUBROUTINE BUILD_IPSTATUS
      USE DATAPOOL
      IMPLICIT NONE
      integer IP
      allocate(IPstatus(MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('BUILD_IPSTATUS, step 1')
#ifndef MPI_PARALL_GRID
      IPstatus=1
#else
      IPstatus=0
      DO IP=1,NP_RES
        IF(ASSOCIATED(IPGL(IPLG(IP))%NEXT)) THEN !interface nodes
          IF(IPGL(IPLG(ip))%NEXT%RANK .ge. MYRANK) THEN  ! interface node is not in the sum already ...
            IPstatus(IP)=1
          ENDIF 
        ELSE
          IPstatus(IP)=1
        ENDIF
      END DO
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE BUILD_TRIANGLE_CORRESPONDENCES
      USE DATAPOOL
      IMPLICIT NONE
      integer :: ListMapped(ne_total)
      integer :: ListMappedB(ne_total)
#ifdef MPI_PARALL_GRID
      integer :: ListFirst(nproc)
      integer :: ListCommon_recv(nproc)
      integer :: ListCommon_send(nproc)
#endif
      integer :: IEfound, IE2, IE, J, IP
      integer :: IPglob, I, MAXMNECCON_TOTAL, idx
      integer :: iProc, MNE_loc, MNEextent_loc
      integer :: nbCommon_send, nbCommon_recv
      integer :: idx_send, idx_recv, iNeigh
      integer :: nbCommon, IPloc, IE_glob, IEloc
      integer :: IP1, IP2, IP3
      integer :: IPglob1, IPglob2, IPglob3
      integer :: iRank, sumExtent, IEadj, eVal
      integer :: I1, I2, J1, J2
      integer, allocatable :: INDX_IE(:,:)
      integer, allocatable :: IE_LocalGlobal(:), StatusNeed(:)
      integer, allocatable :: CCON_total(:)
      integer, allocatable :: eInt(:), dspl_send(:), dspl_recv(:)
      integer, allocatable :: INDXextent_IE(:)
      integer, allocatable :: IEneighbor_V1(:,:)
      integer, allocatable :: IEmembership(:)
      allocate(CCON_total(np_total), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('BUILD_TRIANGLE error 1')
      CCON_TOTAL=0 !global nnegb
      DO IE=1,NE_TOTAL
        DO I=1,3
          IPglob=INEtotal(I,IE)
          CCON_TOTAL(IPglob)=CCON_TOTAL(IPglob) + 1
        END DO
      END DO
      MAXMNECCON_TOTAL=maxval(CCON_TOTAL)
      !
      allocate(INDX_IE(MAXMNECCON_TOTAL, np_total), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('BUILD_TRIANGLE error 2')
      CCON_TOTAL=0
      DO IE=1,NE_TOTAL
        DO I=1,3
          IPglob=INEtotal(I,IE)
          idx=CCON_TOTAL(IPglob)
          CCON_TOTAL(IPglob)=idx + 1
          INDX_IE(idx+1, IPglob)=IE !global ine_gb
        END DO
      END DO
      !
      allocate(IE_LocalGlobal(MNE), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('BUILD_TRIANGLE error 3')
      DO IE=1,MNE
        IP1=INE(1,IE)
        IP2=INE(2,IE)
        IP3=INE(3,IE)
#ifdef MPI_PARALL_GRID
        IPglob1=iplg(IP1)
        IPglob2=iplg(IP2)
        IPglob3=iplg(IP3)
#else
        IPglob1=IP1
        IPglob2=IP2
        IPglob3=IP3
#endif
        IEfound=-1
        DO J=1,CCON_TOTAL(IPglob1)
          IE2=INDX_IE(J, IPglob1)
          IF ((IPglob1 .eq. INEtotal(1,IE2)).and.(IPglob2 .eq. INEtotal(2,IE2)).and.(IPglob3 .eq. INEtotal(3,IE2))) THEN
            IEfound=IE2
          END IF
        END DO
        IF (IEfound .eq. -1) THEN
          write(DBG%FHNDL,*) 'IE=', IE
          write(DBG%FHNDL,*) 'IP123=', IP1, IP2, IP3
          write(DBG%FHNDL,*) 'IPglob123=', IPglob1, IPglob2, IPglob3
          write(DBG%FHNDL,*) 'CCON_TOTAL=', CCON_TOTAL(IPglob1)
          DO J=1,CCON_TOTAL(IPglob1)
            IE2=INDX_IE(J, IPglob1)
            write(DBG%FHNDL,*) 'J=', J, 'IE2=', IE2
            write(DBG%FHNDL,*) 'INE(:,IE2)=', INEtotal(1,IE2), INEtotal(2,IE2), INEtotal(3,IE2)
          END DO
          CALL WWM_ABORT('Did not find the triangle')
        END IF
        IE_LocalGlobal(IE)=IEfound
      END DO
      !
      allocate(IEneighbor_V1(3,MNE), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('BUILD_TRIANGLE error 4')
      DO IE=1,MNE
        DO I1=1,3
          IF (I1 .eq. 3) THEN
            I2=1
          ELSE
            I2=I1+1
          END IF
          IP1=INE(I1,IE)
          IP2=INE(I2,IE)
#ifdef MPI_PARALL_GRID
          IPglob1=iplg(IP1)
          IPglob2=iplg(IP2)
#else
          IPglob1=IP1
          IPglob2=IP2
#endif
          IEadj=-1
          DO J=1,CCON_TOTAL(IPglob1)
            IF (IEadj .eq. -1) THEN
              IE2=INDX_IE(J, IPglob1)
              DO J1=1,3
                IF (J1 .eq. 3) THEN
                  J2=1
                ELSE
                  J2=J1+1
                END IF
                IF ((INEtotal(J1,IE2) .eq. IPglob2).and.(INEtotal(J2,IE2) .eq. IPglob1)) THEN
                  IEadj=IE2
                END IF
              END DO
            END IF            
          END DO
          IEneighbor_V1(I1,IE)=IEadj
        END DO
      END DO
      !
      allocate(StatusNeed(ne_total), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('BUILD_TRIANGLE error 5')
      StatusNeed=0
      DO IE=1,MNE
        IE_glob=IE_LocalGlobal(IE)
        StatusNeed(IE_glob)=1
        DO I=1,3
          IEadj=IEneighbor_V1(I,IE)
          IF (IEadj .ne. -1) THEN
            StatusNeed(IEadj)=1
          END IF
        END DO
      END DO
      MNEextent=sum(StatusNeed)
!      WRITE(STAT%FHNDL,*) 'MNE/MNEextent=', MNE, MNEextent
!      FLUSH(STAT%FHNDL)
      allocate(INDXextent_IE(MNEextent), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('BUILD_TRIANGLE error 6')
      DO IE=1,MNE
        IE_glob=IE_LocalGlobal(IE)
        INDXextent_IE(IE)=IE_glob
        StatusNeed(IE_glob)=0
      END DO
      idx=MNE
      DO IE=1,NE_TOTAL
        IF (StatusNeed(IE) .eq. 1) THEN
          idx=idx+1
          INDXextent_IE(idx)=IE
        END IF
      END DO
      StatusNeed=0
      DO IE=1,MNEextent
        IE_glob=INDXextent_IE(IE)
        StatusNeed(IE_glob)=IE
      END DO
      !
      ! Now building the neighbor list
      !
      allocate(IEneighbor(3,MNE), stat=istat)
      DO IE=1,MNE
        DO I=1,3
          IE_glob=IEneighbor_V1(I,IE)
          IF (IE_glob .eq. -1) THEN
            IEloc=-1
          ELSE
            IEloc=StatusNeed(IE_glob)
          END IF
          IEneighbor(I,IE)=IEloc
        END DO
      END DO
#ifdef MPI_PARALL_GRID
      !
      ! Collecting the ListMNE, ListMNEextent
      !
      allocate(ListMNE(nproc), ListMNEextent(nproc), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('BUILD_TRIANGLE error 7')
      allocate(eInt(2), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('BUILD_TRIANGLE error 8')
      IF (myrank .eq. 0) THEN
        ListMNE(1)=MNE
        ListMNEextent(1)=MNEextent
        DO IPROC=2,nproc
          iRank=IPROC-1
          CALL MPI_RECV(eInt, 2, itype,iRank,2049,comm,istatus,ierr)
          ListMNE(IPROC)=eInt(1)
          ListMNEextent(IPROC)=eInt(2)
        END DO
        DO IPROC=2,nproc
          iRank=IPROC-1
          CALL MPI_SEND(ListMNE,nproc,itype,iRank,2050,comm,ierr)
          CALL MPI_SEND(ListMNEextent,nproc,itype,iRank,2051,comm,ierr)
        END DO
      ELSE
        eInt(1)=MNE
        eInt(2)=MNEextent
        CALL MPI_SEND(eInt,2,itype,0,2049,comm,ierr)
        CALL MPI_RECV(ListMNE, nproc, itype,0,2050,comm,istatus,ierr)
        CALL MPI_RECV(ListMNEextent, nproc, itype,0,2051,comm,istatus,ierr)
      END IF
      deallocate(eInt)
      !
      ! Collecting the ListINDXextent_IE
      !
      sumExtent=sum(ListMNEextent)
      allocate(ListINDXextent_IE(sumExtent), stat=istat)
      IF (myrank .eq. 0) THEN
        idx=0
        DO J=1,MNEextent
          idx=idx+1
          ListINDXextent_IE(idx)=INDXextent_IE(J)
        END DO
        DO IPROC=2,nproc
          iRank=IPROC-1
          MNEextent_loc=ListMNEextent(IPROC)
          allocate(eInt(MNEextent_loc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('BUILD_TRIANGLE error 9')
          CALL MPI_RECV(eInt, MNEextent_loc, itype,iRank,2051,comm,istatus,ierr)
          DO J=1,MNEextent_loc
            idx=idx+1
            ListINDXextent_IE(idx)=eInt(J)
          END DO
          deallocate(eInt)
        END DO
        DO IPROC=2,nproc
          iRank=IPROC-1
          CALL MPI_SEND(ListINDXextent_IE,sumExtent,itype,iRank,2052,comm,ierr)
        END DO
      ELSE
        CALL MPI_SEND(INDXextent_IE,MNEextent,itype,0,2051,comm,ierr)
        CALL MPI_RECV(ListINDXextent_IE, sumExtent, itype,0,2052,comm,istatus,ierr)
      END IF
      !
      ! Now building ListFirst
      !
      ListFirst=0
      DO iProc=2,nproc
        ListFirst(iProc)=ListFirst(iProc-1) + ListMNEextent(iProc-1)
      END DO
#endif
      !
      ! Determine IE membership
      !
      allocate(IEstatus(MNE))
#ifdef MPI_PARALL_GRID
      allocate(IEmembership(ne_total), stat=istat)
      IEmembership=0
      DO iProc=1,nproc
        MNE_loc=ListMNE(iProc)
        DO IE=1,MNE_loc
          IE_glob=ListINDXextent_IE(IE+ListFirst(iProc))
          IF (IEmembership(IE_glob) .eq. 0) THEN
            IEmembership(IE_glob)=iProc
          END IF
        END DO
      END DO
      DO IE=1,MNE
        IE_glob=INDXextent_IE(IE)
        eVal=0
        IF (IEmembership(IE_glob) .eq. myrank+1) eVal=1
        IEstatus(IE)=eVal
      END DO
      deallocate(IEmembership)
#else
      IEstatus=1      
#endif
      WRITE(STAT%FHNDL,*) 'sum(IEstatus)=', sum(IEstatus)
#ifdef MPI_PARALL_GRID
      !
      ! Now building synchronization arrays
      !
      ListMapped=0
      DO IE=1,MNEextent
        IE_glob=INDXextent_IE(IE)
        ListMapped(IE_glob)=IE
      END DO
      ListCommon_send=0
      ListCommon_recv=0
      ie_nnbr_send=0
      ie_nnbr_recv=0
      DO iProc=1,nproc
        IF (iProc .ne. myrank+1) THEN
          MNE_loc=ListMNE(iProc)
          MNEextent_loc=ListMNEextent(iProc)
          ListMappedB=0
          DO IE=1,MNEextent_loc
            IE_glob=ListINDXextent_IE(IE+ListFirst(iProc))
            ListMappedB(IE_glob)=IE
          END DO
          !
          nbCommon_recv=0
          DO IE=1,MNE_loc
            IE_glob=ListINDXextent_IE(IE+ListFirst(iProc))
            IF (ListMapped(IE_glob).gt.0) THEN
              nbCommon_recv=nbCommon_recv+1
            END IF
          END DO
          IF (nbCommon_recv .gt. 0) THEN
            ie_nnbr_recv=ie_nnbr_recv+1
            ListCommon_recv(iProc)=nbCommon_recv
          END IF
          !
          nbCommon_send=0
          DO IE=1,MNE
            IE_glob=INDXextent_IE(IE)
            IF (ListMappedB(IE_glob).gt.0) THEN
              nbCommon_send=nbCommon_send+1
            END IF
          END DO
          IF (nbCommon_send .gt. 0) THEN
            ie_nnbr_send=ie_nnbr_send+1
            ListCommon_send(iProc)=nbCommon_send
          END IF
        END IF
      END DO
      !
      ! Building list of neighbors
      !
      allocate(ListNeigh_ie_send(ie_nnbr_send), ListNeigh_ie_recv(ie_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 27')
      idx_send=0
      idx_recv=0
      DO iProc=1,nproc
        IF (ListCommon_send(iProc) .gt. 0) THEN
          idx_send=idx_send+1
          ListNeigh_ie_send(idx_send)=iProc-1
        END IF
        IF (ListCommon_recv(iProc) .gt. 0) THEN
          idx_recv=idx_recv+1
          ListNeigh_ie_recv(idx_recv)=iProc-1
        END IF
      END DO
      !
      ! Building MPI arrays
      !
      allocate(ie_send_rqst(ie_nnbr_send), ie_recv_rqst(ie_nnbr_recv), ie_send_stat(MPI_STATUS_SIZE,ie_nnbr_send), ie_recv_stat(MPI_STATUS_SIZE,ie_nnbr_recv), ie_send_type(ie_nnbr_send), ie_recv_type(ie_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 38')
      DO iNeigh=1,ie_nnbr_send
        iProc=ListNeigh_ie_send(iNeigh)+1
        nbCommon=ListCommon_send(iProc)
        MNEextent_loc=ListMNEextent(iProc)
        ListMappedB=0
        DO IE=1,MNEextent_loc
          IE_glob=ListINDXextent_IE(IE+ListFirst(iProc))
          ListMappedB(IE_glob)=IE
        END DO
        allocate(dspl_send(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 40')
        idx=0
        DO IE=1,MNE
          IE_glob=INDXextent_IE(IE)
          IEloc=ListMappedB(IE_glob)
          IF (IEloc.gt.0) THEN
            idx=idx+1
            dspl_send(idx)=IE-1
          END IF
        END DO
        call mpi_type_create_indexed_block(nbCommon,1,dspl_send,rtype,ie_send_type(iNeigh), ierr)
        call mpi_type_commit(ie_send_type(iNeigh), ierr)
        deallocate(dspl_send)
      END DO
      DO iNeigh=1,ie_nnbr_recv
        iProc=ListNeigh_ie_recv(iNeigh)+1
        MNE_loc=ListMNE(iProc)
        MNEextent_loc=ListMNEextent(iProc)
        nbCommon=ListCommon_recv(iProc)
        allocate(dspl_recv(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 42')
        idx=0
        DO IE=1,MNE_loc
          IE_glob=ListINDXextent_IE(IE+ListFirst(iProc))
          IEloc=ListMapped(IE_glob)
          IF (IEloc .gt. 0) THEN
            idx=idx+1
            dspl_recv(idx)=IEloc-1
          END IF
        END DO
        call mpi_type_create_indexed_block(nbCommon,1,dspl_recv,rtype,ie_recv_type(iNeigh), ierr)
        call mpi_type_commit(ie_recv_type(iNeigh), ierr)
        deallocate(dspl_recv)
      END DO
      deallocate(INDX_IE, IE_LocalGlobal, StatusNeed, CCON_total, INDXextent_IE, IEneighbor_V1)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE TOTAL_SUMMATION_SCALAR(F, eSum)
      USE DATAPOOL
      IMPLICIT NONE
      real(rkind), intent(in) :: F(MNP)
      real(rkind), intent(out) :: eSum
      real(rkind) eField(1), Lsum(1)
      integer IP, iProc
      eSum=ZERO
      DO IP=1,MNP
        IF (IPstatus(IP) .eq. 1) THEN
          eSum = eSum + F(IP)
        END IF
      END DO
!      WRITE(STAT%FHNDL,*) 'sum(AC)=', sum(AC)
!      WRITE(STAT%FHNDL,*) 'After summation, sum(Lsum)=', sum(Lsum)
#ifdef MPI_PARALL_GRID
      IF (myrank == 0) THEN
        DO iProc=2,nproc
          CALL MPI_RECV(eField,1,rtype, iProc-1, 43, comm, istatus, ierr)
          eSum=eSum + eField(1)
        END DO
        Lsum(1)=eSum
        DO iProc=2,nproc
          CALL MPI_SEND(Lsum,1,rtype, iProc-1, 13, comm, ierr)
        END DO
      ELSE
        Lsum(1)=eSum
        CALL MPI_SEND(Lsum,1,rtype, 0, 43, comm, ierr)
        CALL MPI_RECV(Lsum,1,rtype, 0, 13, comm, istatus, ierr)
        eSum=Lsum(1)
      END IF
#endif
!      WRITE(STAT%FHNDL,*) 'At leaving, sum(Lsum)=', sum(Lsum)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE TOTAL_SUMMATION_AC(AC, Lsum)
      USE DATAPOOL
      IMPLICIT NONE
      real(rkind), intent(in) :: AC(MSC, MDC, MNP)
      real(rkind), intent(out) :: Lsum(MSC, MDC)
      real(rkind) eField(MSC, MDC)
      integer IP, iProc
      Lsum=ZERO
      DO IP=1,MNP
        IF (IPstatus(IP) .eq. 1) THEN
          Lsum = Lsum + AC(:,:,IP)
        END IF
      END DO
!      WRITE(STAT%FHNDL,*) 'sum(AC)=', sum(AC)
!      WRITE(STAT%FHNDL,*) 'After summation, sum(Lsum)=', sum(Lsum)
#ifdef MPI_PARALL_GRID
      IF (myrank == 0) THEN
        DO iProc=2,nproc
          CALL MPI_RECV(eField,MSC*MDC,rtype, iProc-1, 43, comm, istatus, ierr)
          Lsum=Lsum + eField
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(Lsum,MSC*MDC,rtype, iProc-1, 13, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(Lsum,MSC*MDC,rtype, 0, 43, comm, ierr)
        CALL MPI_RECV(Lsum,MSC*MDC,rtype, 0, 13, comm, istatus, ierr)
      END IF
#endif     
!      WRITE(STAT%FHNDL,*) 'At leaving, sum(Lsum)=', sum(Lsum)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE Print_SumAC2(string)
      USE DATAPOOL
      implicit NONE
      character(len=*), intent(in) :: string
      real(rkind) :: Lsum(MSC,MDC)
!      WRITE(STAT%FHNDL,*) 'Direct sum(AC2)=', sum(AC2)
      CALL TOTAL_SUMMATION_AC(AC2, Lsum)
      WRITE(STAT%FHNDL,*) 'sum(AC2)=', sum(Lsum),' at step:', TRIM(string)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE Print_SumScalar(F, string)
      USE DATAPOOL
      implicit NONE
      real(rkind), intent(in) :: F(MNP)
      character(len=*) :: string
      real(rkind) :: eSum
      CALL TOTAL_SUMMATION_SCALAR(F, eSum)
      WRITE(STAT%FHNDL,*) 'sum(F)=', eSum,' mesg:', TRIM(string)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_TOTAL_COHERENCY_ERROR(MSCeffect, ACw, Lerror)
      USE DATAPOOL
      implicit none
      integer, intent(in) :: MSCeffect
      real(rkind), intent(in) :: ACw(MSCeffect, MDC, MNP)
      real(rkind), intent(out) :: Lerror
      real(rkind), allocatable :: ACtotal(:,:,:)
      real(rkind), allocatable :: ACloc(:,:,:)
      real(rkind) :: rbuf_real(1)
      integer, allocatable :: ListFirstMNP(:)
      integer, allocatable :: eStatus(:)
      integer IP, iProc, IPglob, IS, ID
      integer MNPloc
#ifdef MPI_PARALL_GRID
      IF (myrank == 0) THEN
        Lerror=0
        allocate(ListFirstMNP(nproc), eStatus(np_global), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 69')
        ListFirstMNP=0
        eStatus=0
        DO iProc=2,nproc
          ListFirstMNP(iProc)=ListFirstMNP(iProc-1) + ListMNP(iProc-1)
        END DO
        allocate(ACtotal(MSCeffect, MDC, np_global), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 70')
        DO IP=1,MNP
          IPglob=iplg(IP)
          ACtotal(:,:,IPglob)=ACw(:,:,IP)
          eStatus(IPglob)=1
        END DO
        DO iProc=2,nproc
          MNPloc=ListMNP(iProc)
          allocate(ACloc(MSCeffect, MDC, MNPloc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 71')
          CALL MPI_RECV(ACloc,MNPloc*MSCeffect*MDC,rtype, iProc-1, 53, comm, istatus, ierr)
          DO IP=1,MNPloc
            IPglob=ListIPLG(IP+ListFirstMNP(iProc))
            IF (eStatus(IPglob) == 1) THEN
              DO IS=1,MSCeffect
                DO ID=1,MDC
                  Lerror=Lerror+abs(ACtotal(IS,ID,IPglob)-ACloc(IS,ID,IP))
                END DO
              END DO
            ELSE
              eStatus(IPglob)=1
              ACtotal(:,:,IPglob)=ACloc(:,:,IP)
            END IF
          END DO
          deallocate(ACloc)
        END DO
        deallocate(ListFirstMNP, ACtotal, eStatus)
        rbuf_real(1)=Lerror
        DO iProc=2,nproc
          CALL MPI_SEND(rbuf_real,1,rtype, iProc-1, 23, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(ACw,MNP*MSCeffect*MDC,rtype, 0, 53, comm, ierr)
        CALL MPI_RECV(rbuf_real,1,rtype, 0, 23, comm, istatus, ierr)
        Lerror=rbuf_real(1)
      END IF
#else
      Lerror=0 ! in serial the coherency error is zero
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_TOTAL_COHERENCY_ERROR_NPRES(MSCeffect, ACw, Lerror)
      USE DATAPOOL
      implicit none
      integer, intent(in) :: MSCeffect
      real(rkind), intent(in) :: ACw(MSCeffect, MDC, MNP)
      real(rkind), intent(out) :: Lerror
      real(rkind), allocatable :: ACtotal(:,:,:)
      real(rkind), allocatable :: ACloc(:,:,:)
      real(rkind) :: rbuf_real(1)
      integer, allocatable :: ListFirstMNP(:)
      integer, allocatable :: eStatus(:)
      integer IP, iProc, IPglob, IS, ID
      integer NP_RESloc
#ifdef MPI_PARALL_GRID
      IF (myrank == 0) THEN
        Lerror=0
        allocate(ListFirstMNP(nproc), eStatus(np_global), ACtotal(MSCeffect, MDC, np_global), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 72')
        ListFirstMNP=0
        eStatus=0
        DO iProc=2,nproc
          ListFirstMNP(iProc)=ListFirstMNP(iProc-1) + ListMNP(iProc-1)
        END DO
        DO IP=1,NP_RES
          IPglob=iplg(IP)
          ACtotal(:,:,IPglob)=ACw(:,:,IP)
          eStatus(IPglob)=1
        END DO
        DO iProc=2,nproc
          NP_RESloc=ListNP_RES(iProc)
          allocate(ACloc(MSCeffect, MDC, NP_RESloc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 73')
          CALL MPI_RECV(ACloc,MSCeffect*MDC*NP_RESloc,rtype, iProc-1, 53, comm, istatus, ierr)
          DO IP=1,NP_RESloc
            IPglob=ListIPLG(IP+ListFirstMNP(iProc))
            IF (eStatus(IPglob) == 1) THEN
              DO IS=1,MSCeffect
                DO ID=1,MDC
                  Lerror=Lerror+abs(ACtotal(IS,ID,IPglob)-ACloc(IS,ID,IP))
                END DO
              END DO
            ELSE
              eStatus(IPglob)=1
              ACtotal(:,:,IPglob)=ACloc(:,:,IP)
            END IF
          END DO
          deallocate(ACloc)
        END DO
        deallocate(ListFirstMNP, ACtotal, eStatus)
        rbuf_real(1)=Lerror
        DO iProc=2,nproc
          CALL MPI_SEND(rbuf_real,1,rtype, iProc-1, 23, comm, ierr)
        END DO
      ELSE
        allocate(ACloc(MSCeffect, MDC, NP_RES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 74')
        DO IP=1,NP_RES
          ACloc(:,:,IP)=ACw(:,:,IP)
        END DO
        CALL MPI_SEND(ACloc,NP_RES*MSCeffect*MDC,rtype, 0, 53, comm, ierr)
        deallocate(ACloc)
        CALL MPI_RECV(rbuf_real,1,rtype, 0, 23, comm, istatus, ierr)
        Lerror=rbuf_real(1)
      END IF
#else
      Lerror=0 ! in serial the coherency error is zero
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef MPI_PARALL_GRID
      SUBROUTINE COLLECT_ALL_IPLG
      USE DATAPOOL
      implicit none
      integer, allocatable :: rbuf_int(:)
      integer len, iProc, IP, idx, sumMNP
      allocate(ListMNP(nproc), ListNP_RES(nproc), rbuf_int(2), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 15')
      IF (myrank == 0) THEN
        ListMNP(1)=MNP
        ListNP_RES(1)=NP_RES
        DO iProc=2,nproc
          CALL MPI_RECV(rbuf_int,2,itype, iProc-1, 257, comm, istatus, ierr)
          ListMNP(iProc)=rbuf_int(1)
          ListNP_RES(iProc)=rbuf_int(2)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(ListMNP,nproc,itype, iProc-1, 263, comm, ierr)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(ListNP_RES,nproc,itype, iProc-1, 571, comm, ierr)
        END DO
      ELSE
        rbuf_int(1)=MNP
        rbuf_int(2)=NP_RES
        CALL MPI_SEND(rbuf_int,2,itype, 0, 257, comm, ierr)
        CALL MPI_RECV(ListMNP,nproc,itype, 0, 263, comm, istatus, ierr)
        CALL MPI_RECV(ListNP_RES,nproc,itype, 0, 571, comm, istatus, ierr)
      END IF
      deallocate(rbuf_int)
      sumMNP=sum(ListMNP)
      allocate(ListIPLG(sumMNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 16')
      IF (myrank == 0) THEN
        idx=0
        DO IP=1,MNP
          idx=idx+1
          ListIPLG(idx)=iplg(IP)
        END DO
        DO iProc=2,nproc
          len=ListMNP(iProc)
          allocate(rbuf_int(len), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 17')
          CALL MPI_RECV(rbuf_int,len,itype, iProc-1, 269, comm, istatus, ierr)
          DO IP=1,len
            idx=idx+1
            ListIPLG(idx)=rbuf_int(IP)
          END DO
          deallocate(rbuf_int)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(ListIPLG,sumMNP,itype, iProc-1, 271, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(iplg,MNP,itype, 0, 269, comm, ierr)
        CALL MPI_RECV(ListIPLG,sumMNP,itype, 0, 271, comm, istatus, ierr)
      END IF
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef MPI_PARALL_GRID
      SUBROUTINE COLLECT_ALL_IA_JA
      USE DATAPOOL
      implicit none
      integer, allocatable :: rbuf_int(:)
      integer len, iProc, IP, idx
      integer sumIAsiz, sumNNZ
      allocate(ListNNZ(nproc), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 20')
      !
      ! Collecting NNZ
      !
      IF (myrank == 0) THEN
        ListNNZ(1)=NNZ
        allocate(rbuf_int(1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 21')
        DO iProc=2,nproc
          CALL MPI_RECV(rbuf_int,1,itype, iProc-1, 257, comm, istatus, ierr)
          ListNNZ(iProc)=rbuf_int(1)
        END DO
        deallocate(rbuf_int)
        DO iProc=2,nproc
          CALL MPI_SEND(ListNNZ,nproc,itype, iProc-1, 263, comm, ierr)
        END DO
      ELSE
        allocate(rbuf_int(1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 22')
        rbuf_int(1)=NNZ
        CALL MPI_SEND(rbuf_int,1,itype, 0, 257, comm, ierr)
        deallocate(rbuf_int)
        CALL MPI_RECV(ListNNZ,nproc,itype, 0, 263, comm, istatus, ierr)
      END IF
      !
      ! Collecting IA
      !
      sumIAsiz=sum(ListMNP) + nproc
      allocate(ListIA(sumIAsiz), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 23')
      IF (myrank == 0) THEN
        idx=0
        DO IP=1,MNP+1
          idx=idx+1
          ListIA(idx)=IA(IP)
        END DO
        DO iProc=2,nproc
          len=ListMNP(iProc)+1
          allocate(rbuf_int(len), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 24')
          CALL MPI_RECV(rbuf_int,len,itype, iProc-1, 269, comm, istatus, ierr)
          DO IP=1,len
            idx=idx+1
            ListIA(idx)=rbuf_int(IP)
          END DO
          deallocate(rbuf_int)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(ListIA,sumIAsiz,itype, iProc-1, 271, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(IA,MNP+1,itype, 0, 269, comm, ierr)
        CALL MPI_RECV(ListIA,sumIAsiz,itype, 0, 271, comm, istatus, ierr)
      END IF
      !
      ! Collecting JA
      !
      sumNNZ=sum(ListNNZ)
      allocate(ListJA(sumNNZ), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 25')
      IF (myrank == 0) THEN
        idx=0
        DO IP=1,NNZ
          idx=idx+1
          ListJA(idx)=JA(IP)
        END DO
        DO iProc=2,nproc
          len=ListNNZ(iProc)
          allocate(rbuf_int(len), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 26')
          CALL MPI_RECV(rbuf_int,len,itype, iProc-1, 569, comm, istatus, ierr)
          DO IP=1,len
            idx=idx+1
            ListJA(idx)=rbuf_int(IP)
          END DO
          deallocate(rbuf_int)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(ListJA,sumNNZ,itype, iProc-1, 467, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(JA,NNZ,itype, 0, 569, comm, ierr)
        CALL MPI_RECV(ListJA,sumNNZ,itype, 0, 467, comm, istatus, ierr)
      END IF
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef MPI_PARALL_GRID
      SUBROUTINE SYMM_GRAPH_BUILD_ADJ(AdjGraph)
      USE DATAPOOL
      implicit none
      type(Graph), intent(inout) :: AdjGraph
      CALL KERNEL_GRAPH_BUILD_ADJ(AdjGraph, wwm_nnbr, wwm_ListNeigh)
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef MPI_PARALL_GRID
      SUBROUTINE GRAPH_BUILD_ADJ(AdjGraph)
      USE datapool, only : nnbr_p, nbrrank_p, Graph
      implicit none
      type(Graph), intent(inout) :: AdjGraph
      integer :: ListNe(nnbr_p)
      integer I
      DO I=1,nnbr_p
        ListNe(I)=nbrrank_p(I)+1
      END DO
      CALL KERNEL_GRAPH_BUILD_ADJ(AdjGraph, nnbr_p, ListNe)
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef MPI_PARALL_GRID
      SUBROUTINE KERNEL_GRAPH_BUILD_ADJ(AdjGraph, nb, ListNe)
      USE DATAPOOL
      implicit none
      integer, intent(in) :: nb
      integer, intent(in) :: ListNe(nb)
      type(Graph), intent(inout) :: AdjGraph
      integer, allocatable :: rbuf_int(:)
      integer I, iProc
      integer idx, eDeg, nbEdge, iEdge
      AdjGraph % nbVert=nproc
      IF (myrank.eq.0) THEN
        allocate(AdjGraph % ListDegree(nproc), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 1')
        AdjGraph % ListDegree(1)=nb
        allocate(rbuf_int(1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 2')
        DO iProc=2,nproc
          CALL MPI_RECV(rbuf_int,1,itype, iProc-1, 19, comm, istatus, ierr)
          AdjGraph % ListDegree(iProc)=rbuf_int(1)
        END DO
        deallocate(rbuf_int)
        nbEdge=0
        DO iProc=1,nproc
          nbEdge=nbEdge + AdjGraph % ListDegree(iProc)
        END DO
        AdjGraph % nbEdge=nbEdge
        allocate(AdjGraph % ListEdge(nbEdge,2), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 3')
        idx=0
        eDeg=AdjGraph % ListDegree(1)
        DO I=1,eDeg
          idx=idx+1
          AdjGraph % ListEdge(idx,1)=1
          AdjGraph % ListEdge(idx,2)=ListNe(I)
        END DO
        DO iProc=2,nproc
          eDeg=AdjGraph % ListDegree(iProc)
          allocate(rbuf_int(eDeg), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 4')
          CALL MPI_RECV(rbuf_int,eDeg,itype, iProc-1, 24, comm, istatus, ierr)
          DO I=1,eDeg
            idx=idx+1
            AdjGraph % ListEdge(idx,1)=iProc
            AdjGraph % ListEdge(idx,2)=rbuf_int(I)
          END DO
          deallocate(rbuf_int)
        END DO
        allocate(rbuf_int(1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 5')
        rbuf_int(1)=nbEdge
        DO iProc=2,nproc
          CALL MPI_SEND(rbuf_int,1,itype, iProc-1, 30, comm, ierr)
        END DO
        deallocate(rbuf_int)
        !
        allocate(rbuf_int(nproc), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 6')
        DO iProc=1,nproc
          rbuf_int(iProc)=AdjGraph % ListDegree(iProc)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(rbuf_int,nproc,itype, iProc-1, 32, comm, ierr)
        END DO
        deallocate(rbuf_int)
        !
        allocate(rbuf_int(nbEdge*2), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 7')
        DO iEdge=1,nbEdge
          rbuf_int(2*(iEdge-1)+1)=AdjGraph % ListEdge(iEdge,1)
          rbuf_int(2*(iEdge-1)+2)=AdjGraph % ListEdge(iEdge,2)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(rbuf_int,2*nbEdge,itype, iProc-1, 34, comm, ierr)
        END DO
      ELSE
        allocate(rbuf_int(1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 8')
        rbuf_int(1)=nb
        CALL MPI_SEND(rbuf_int,1,itype, 0, 19, comm, ierr)
        deallocate(rbuf_int)
        !
        allocate(rbuf_int(nb), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 9')
        DO I=1,nb
          rbuf_int(I)=ListNe(I)
        END DO
        CALL MPI_SEND(rbuf_int,nb,itype, 0, 24, comm, ierr)
        deallocate(rbuf_int)
        !
        allocate(rbuf_int(1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 10')
        CALL MPI_RECV(rbuf_int,1,itype, 0, 30, comm, istatus, ierr)
        nbEdge=rbuf_int(1)
        deallocate(rbuf_int)
        AdjGraph % nbEdge=nbEdge
        !
        allocate(rbuf_int(nproc), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 11')
        CALL MPI_RECV(rbuf_int,nproc,itype, 0, 32, comm, istatus, ierr)
        allocate(AdjGraph % ListDegree(nproc), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 12')
        AdjGraph % ListDegree=rbuf_int
        deallocate(rbuf_int)
        !
        allocate(rbuf_int(2*nbEdge), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 13')
        CALL MPI_RECV(rbuf_int,2*nbEdge,itype, 0, 34, comm, istatus, ierr)
        allocate(AdjGraph % ListEdge(nbEdge,2), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 14')
        DO iEdge=1,nbEdge
          AdjGraph % ListEdge(iEdge,1)=rbuf_int(2*(iEdge-1)+1)
          AdjGraph % ListEdge(iEdge,2)=rbuf_int(2*(iEdge-1)+2)
        END DO
        deallocate(rbuf_int)
      ENDIF
      AdjGraph % MaxDeg=maxval(AdjGraph % ListDegree)
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GRAPH_TEST_CONNECT(AdjGraph, result)
      USE DATAPOOL
      implicit none
      type(Graph), intent(in) :: AdjGraph
      integer, intent(out) :: result
      integer :: ListStatus(AdjGraph%nbVert)
      integer :: ListPosFirst(AdjGraph%nbVert)
      integer idx, iVert, nbVert, eAdj
      integer eDeg, sizConn, I, IsFinished
      idx=0
      nbVert=AdjGraph%nbVert
      ListStatus=0
      DO iVert=1,nbVert
        ListPosFirst(iVert)=idx
        idx=idx+AdjGraph % ListDegree(iVert)
      END DO
      ListStatus(1)=1
      DO
        IsFinished=1
        DO iVert=1,nbVert
          IF (ListStatus(iVert) == 1) THEN
            eDeg=AdjGraph%ListDegree(iVert)
            idx=ListPosFirst(iVert)
            DO I=1,eDeg
              eAdj=AdjGraph%ListEdge(idx+I,1)
              IF (ListStatus(eAdj) == 0) THEN
                IsFinished=0
              END IF
              ListStatus(eAdj)=1
            END DO
          END IF
        END DO
        IF (IsFinished == 1) THEN
          EXIT
        END IF
      END DO
      sizConn=0
      DO iVert=1,nbVert
        IF (ListStatus(iVert) == 1) THEN
          sizConn=sizConn+1
        END IF
      END DO
      IF (sizConn == nbVert) THEN
        result=1
      END IF
      result=0
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef MPI_PARALL_GRID
      SUBROUTINE SETUP_ONED_SCATTER_ARRAY
      USE DATAPOOL
      IMPLICIT NONE
      integer :: ListFirst(nproc)
      integer MNPloc, iProc, IP, IP_glob
      integer, allocatable :: dspl_oned(:), dspl_twod(:)
      ListFirst=0
      DO iProc=2,nproc
        ListFirst(iProc)=ListFirst(iProc-1) + ListMNP(iProc-1)
      END DO
      IF (myrank .eq. 0) THEN
        allocate(oned_send_rqst(nproc-1), oned_send_stat(MPI_STATUS_SIZE,nproc-1), oned_send_type(nproc-1), stat=istat)
        allocate(twod_send_rqst(nproc-1), twod_send_stat(MPI_STATUS_SIZE,nproc-1), twod_send_type(nproc-1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('error in IOBPtotal allocate')
        DO iProc=2,nproc
          MNPloc=ListMNP(iProc)
          WRITE(STAT%FHNDL,*) 'iProc, MNPloc=', iProc, MNPloc
          allocate(dspl_oned(MNPloc), dspl_twod(MNPloc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('error in IOBPtotal allocate')
          DO IP=1,MNPloc
            IP_glob=ListIPLG(IP+ListFirst(iProc))
            dspl_oned(IP)=IP_glob-1
            dspl_twod(IP)=2*(IP_glob-1)
          END DO
          call mpi_type_create_indexed_block(MNPloc,1,dspl_oned,rtype,oned_send_type(iProc-1), ierr)
          call mpi_type_commit(oned_send_type(iProc-1), ierr)
          call mpi_type_create_indexed_block(MNPloc,2,dspl_twod,rtype,twod_send_type(iProc-1), ierr)
          call mpi_type_commit(twod_send_type(iProc-1), ierr)
          deallocate(dspl_oned, dspl_twod)
        END DO
        FLUSH(STAT%FHNDL)
      END IF
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SCATTER_ONED_ARRAY(Vtotal, Vlocal)
      USE DATAPOOL
      IMPLICIT NONE
      real(rkind), intent(in) :: Vtotal(np_total)
      real(rkind), intent(out) :: Vlocal(MNP)
      integer iProc, IP
#ifdef MPI_PARALL_GRID
      IF (myrank .eq. 0) THEN
        DO iProc=2,nproc
          CALL mpi_isend(Vtotal, 1, oned_send_type(iProc-1), iProc-1, 2030, comm, oned_send_rqst(iProc-1), ierr)
        END DO
        DO IP=1,MNP
          Vlocal(IP)=Vtotal(iplg(IP))
        END DO
        IF (nproc > 1) THEN
          CALL MPI_WAITALL(nproc-1, oned_send_rqst, oned_send_stat, ierr)
        END IF
      ELSE
        CALL MPI_RECV(Vlocal, MNP, rtype, 0, 2030, comm, istatus, ierr)
      END IF
#else
      Vlocal = Vtotal
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SCATTER_TWOD_ARRAY(Vtotal, Vlocal)
      USE DATAPOOL
      IMPLICIT NONE
      real(rkind), intent(in) :: Vtotal(2,np_total)
      real(rkind), intent(out) :: Vlocal(2,MNP)
#ifdef MPI_PARALL_GRID
      integer iProc, IP
      IF (myrank .eq. 0) THEN
        DO iProc=2,nproc
          CALL mpi_isend(Vtotal, 1, twod_send_type(iProc-1), iProc-1, 2068, comm, twod_send_rqst(iProc-1), ierr)
        END DO
        DO IP=1,MNP
          Vlocal(:,IP)=Vtotal(:,iplg(IP))
        END DO
        IF (nproc > 1) THEN
          CALL MPI_WAITALL(nproc-1, twod_send_rqst, twod_send_stat, ierr)
        END IF
      ELSE
        CALL MPI_RECV(Vlocal, 2*MNP, rtype, 0, 2068, comm, istatus, ierr)
      END IF
#else
      Vlocal = Vtotal
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef MPI_PARALL_GRID
      SUBROUTINE SETUP_BOUNDARY_SCATTER_REDUCE_ARRAY
      USE DATAPOOL
      IMPLICIT NONE
      logical, SAVE ::IsSetupBoundaryScatterReduceArrayDone = .FALSE.
      integer :: ListFirst(nproc)
      integer MNPloc, iProc, IP, IP_glob
      integer, allocatable :: dspl_spparm(:), Indexes(:)
      integer, allocatable :: dspl_wbac(:)
      integer :: NbSend(nproc)
      integer irank, eSend, idx, idx_nbproc, eIdx
      IF (IsSetupBoundaryScatterReduceArrayDone .eqv. .TRUE.) THEN
        RETURN
      END IF
      IsSetupBoundaryScatterReduceArrayDone=.TRUE.
      ListFirst=0
      DO iProc=2,nproc
        ListFirst(iProc)=ListFirst(iProc-1) + ListMNP(iProc-1)
      END DO
      bound_nbproc=0
      rank_hasboundary=-1
      DO irank=0,nproc-1
        iProc=irank+1
        eSend=0
        IF (irank .ne. rank_boundary) THEN
          MNPloc=ListMNP(iProc)
          DO IP=1,MNPloc
            IP_glob=ListIPLG(IP+ListFirst(iProc))
            IF ((IOBPtotal(IP_glob) .eq. 2).or.(IOBPtotal(IP_glob) .eq. 4)) THEN
              eSend=eSend+1
            END IF
          END DO
          IF (eSend .gt. 0) THEN
            bound_nbproc=bound_nbproc+1
            rank_hasboundary=irank
          END IF
        END IF
        NbSend(iProc)=eSend
        WRITE(STAT%FHNDL,*) 'iProc=', iProc, ' eSend=', eSend
        FLUSH(STAT%FHNDL)
      END DO
      IF (myrank .eq. rank_boundary) THEN
        WRITE(STAT%FHNDL,*) 'bound_nbproc=', bound_nbproc
        FLUSH(STAT%FHNDL)
        allocate(bound_listproc(bound_nbproc), Indexes(np_total), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('allocate error')
        !
        allocate(spparm_rqst(bound_nbproc), spparm_stat(MPI_STATUS_SIZE,bound_nbproc), spparm_type(bound_nbproc), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('error in IOBPtotal allocate')
        !
        allocate(wbac_rqst(bound_nbproc), wbac_stat(MPI_STATUS_SIZE,bound_nbproc), wbac_type(bound_nbproc), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('error in IOBPtotal allocate')
        !
        idx=0
        DO IP_glob=1,np_total
          IF ((IOBPtotal(IP_glob) .eq. 2).or.(IOBPtotal(IP_glob) .eq. 4)) THEN
            idx=idx+1
            Indexes(IP_glob)=idx
          END IF
        END DO
        WRITE(STAT%FHNDL,*) 'idx=', idx
        FLUSH(STAT%FHNDL)

        IF (IWBMNP .gt. 0) THEN
          allocate(Indexes_boundary(IWBMNP), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('error in IOBPtotal allocate')
          DO IP=1,IWBMNP
            IP_glob=iplg(IWBNDLC(IP))
            Indexes_boundary(IP)=Indexes(IP_glob)
          END DO
        END IF
        idx_nbproc=0
        DO irank=0,nproc-1
          iProc=irank+1
          eSend=NbSend(iProc)
          IF ((irank .ne. rank_boundary).and.(eSend.gt.0)) THEN
            idx_nbproc=idx_nbproc+1
            bound_listproc(idx_nbproc)=irank
            MNPloc=ListMNP(iProc)
            allocate(dspl_spparm(eSend), dspl_wbac(eSend), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('error in IOBPtotal allocate')
            idx=0
            DO IP=1,MNPloc
              IP_glob=ListIPLG(IP+ListFirst(iProc))
              IF ((IOBPtotal(IP_glob) .eq. 2).or.(IOBPtotal(IP_glob) .eq. 4)) THEN
                idx=idx+1
                eIdx=Indexes(IP_glob)
                dspl_spparm(idx)=8*(eIdx-1)
                dspl_wbac(idx)=MSC*MDC*(eIdx-1)
                WRITE(STAT%FHNDL,*) 'idx=', idx, 'eIdx=', eIdx
                FLUSH(STAT%FHNDL)
              END IF
            END DO
            call mpi_type_create_indexed_block(eSend,8,dspl_spparm,rtype,spparm_type(idx_nbproc), ierr)
            call mpi_type_commit(spparm_type(idx_nbproc), ierr)
            call mpi_type_create_indexed_block(eSend,MSC*MDC,dspl_wbac,rtype,wbac_type(idx_nbproc), ierr)
            call mpi_type_commit(wbac_type(idx_nbproc), ierr)
            deallocate(dspl_spparm, dspl_wbac)
          END IF
        END DO
        deallocate(Indexes)
      END IF
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SCATTER_BOUNDARY_ARRAY_SPPARM
      USE DATAPOOL
      IMPLICIT NONE
      integer iProc, IP, irank
      IF ((IWBMNP .eq. 0).and.(myrank.ne.rank_boundary)) THEN
        RETURN
      END IF
#ifdef MPI_PARALL_GRID
      IF (myrank .eq. rank_boundary) THEN
        DO irank=1,bound_nbproc
          CALL mpi_isend(SPPARM_GL, 1, spparm_type(irank), bound_listproc(irank), 2072, comm, spparm_rqst(irank), ierr)
        END DO
        DO IP=1,IWBMNP
          SPPARM(:,IP)=SPPARM_GL(:,Indexes_boundary(IP))
        END DO
        IF (bound_nbproc > 0) THEN
          CALL MPI_WAITALL(bound_nbproc, spparm_rqst, spparm_stat, ierr)
        END IF
      ELSE
        CALL MPI_RECV(SPPARM, 8*IWBMNP, rtype, rank_boundary, 2072, comm, istatus, ierr)
      END IF
#else
      SPPARM = SPPARM_GL
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SCATTER_BOUNDARY_ARRAY_WBAC(WBACOUT)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(OUT)   :: WBACOUT(MSC,MDC,IWBMNP)
      integer iProc, IP, irank
      IF ((IWBMNP .eq. 0).and.(myrank.ne.rank_boundary)) THEN
        RETURN
      END IF
#ifdef MPI_PARALL_GRID
      IF (myrank .eq. rank_boundary) THEN
        DO irank=1,bound_nbproc
          CALL mpi_isend(WBAC_GL, 1, spparm_type(irank), bound_listproc(irank), 2096, comm, spparm_rqst(irank), ierr)
        END DO
        DO IP=1,IWBMNP
          WBACOUT(:,:,IP)=WBAC_GL(:,:,Indexes_boundary(IP))
        END DO
        IF (bound_nbproc > 0) THEN
          CALL MPI_WAITALL(bound_nbproc, spparm_rqst, spparm_stat, ierr)
        END IF
      ELSE
        CALL MPI_RECV(WBACOUT, MSC*MDC*IWBMNP, rtype, rank_boundary, 2096, comm, istatus, ierr)
      END IF
#else
      WBACOUT = WBAC_GL
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE REDUCE_BOUNDARY_ARRAY_SPPARM
      USE DATAPOOL
      IMPLICIT NONE
      integer iProc, IP, irank
#ifndef MPI_PARALL_GRID
      IF (LINHOM) THEN
        SPPARM_GL=SPPARM
      ELSE
        DO IP=1,IWBMNPGL
          SPPARM_GL(:,IP)=SPPARM(:,1)
        END DO
      ENDIF
#else
      IF ((IWBMNP .eq. 0).and.(myrank.ne.rank_boundary)) THEN
        RETURN
      END IF
      IF (LINHOM) THEN
        IF (myrank .eq. rank_boundary) THEN
          DO irank=1,bound_nbproc
            CALL mpi_irecv(SPPARM_GL, 1, spparm_type(irank), bound_listproc(irank), 2099, comm, spparm_rqst(irank), ierr)
          END DO
          DO IP=1,IWBMNP
            SPPARM_GL(:, Indexes_boundary(IP))=SPPARM(:, IP)
          END DO
          IF (bound_nbproc > 0) THEN
            CALL MPI_WAITALL(bound_nbproc, spparm_rqst, spparm_stat, ierr)
          END IF
        ELSE
          CALL MPI_SEND(SPPARM, 8*IWBMNP, rtype, rank_boundary, 2099, comm, ierr)
        END IF
      ELSE
        IF (rank_boundary .ne. rank_hasboundary) THEN
          IF (myrank .eq. rank_hasboundary) THEN
            CALL MPI_SEND(SPPARM,8,rtype, rank_boundary, 2045, comm, ierr)
          END IF
          IF (myrank .eq. rank_boundary) THEN
            CALL MPI_RECV(SPPARM,8,rtype, rank_hasboundary, 2045, comm, istatus, ierr)
          END IF
        END IF
        IF (myrank .eq. rank_boundary) THEN
          DO IP=1,IWBMNPGL
            SPPARM_GL(:,IP)=SPPARM(:,1)
          END DO
        END IF
      END IF
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE REDUCE_BOUNDARY_ARRAY_WBAC
      USE DATAPOOL
      IMPLICIT NONE
      integer iProc, IP, idx_proc
#ifndef MPI_PARALL_GRID
      IF (LINHOM) THEN
        WBAC_GL=WBAC
      ELSE
        DO IP=1,IWBMNPGL
          WBAC_GL(:,:,IP)=WBAC(:,:,1)
        END DO
      ENDIF
#else
      IF ((IWBMNP .eq. 0).and.(myrank.ne.rank_boundary)) THEN
        RETURN
      END IF
      IF (LINHOM) THEN
        IF (myrank .eq. rank_boundary) THEN
          WRITE(STAT%FHNDL,*) 'Before data receiving'
          WRITE(STAT%FHNDL,*) 'IWBMNPGL=', IWBMNPGL
          WRITE(STAT%FHNDL,*) 'MSC/MDC=', MSC,MDC
          WRITE(STAT%FHNDL,*) 'allocated(WBAC_GL)=', allocated(WBAC_GL)
          WRITE(STAT%FHNDL,*) 'size(WBAC_GL)=', size(WBAC_GL)
          FLUSH(STAT%FHNDL)
          WBAC_GL=0
          DO idx_proc=1,bound_nbproc
            WRITE(STAT%FHNDL,*) 'idx_proc/eProc=', idx_proc, bound_listproc(idx_proc)
            FLUSH(STAT%FHNDL)
            CALL mpi_irecv(WBAC_GL, 1, wbac_type(idx_proc), bound_listproc(idx_proc), 2040, comm, wbac_rqst(idx_proc), ierr)
            WRITE(STAT%FHNDL,*) 'MPI_IRECV ierr=', ierr
            FLUSH(STAT%FHNDL)
          END DO
          WRITE(STAT%FHNDL,*) 'IWBMNP=', IWBMNP
          WRITE(STAT%FHNDL,*) 'IWBMNPGL=', IWBMNPGL
          FLUSH(STAT%FHNDL)
          DO IP=1,IWBMNP
            WBAC_GL(:,:,Indexes_boundary(IP))=WBAC(:,:,IP)
          END DO
          IF (bound_nbproc > 0) THEN
            CALL MPI_WAITALL(bound_nbproc, wbac_rqst, wbac_stat, ierr)
            WRITE(STAT%FHNDL,*) 'MPI_WAITALL ierr=', ierr
            FLUSH(STAT%FHNDL)
          END IF
          WRITE(STAT%FHNDL,*) 'sum(WBAC_GL)=', sum(WBAC_GL)
          FLUSH(STAT%FHNDL)
        ELSE
          WRITE(STAT%FHNDL,*) 'Before data sending IWBMNP=', IWBMNP
          WRITE(STAT%FHNDL,*) 'sum(WBAC)=', sum(WBAC)
          FLUSH(STAT%FHNDL)
          CALL MPI_SEND(WBAC, MSC*MDC*IWBMNP, rtype, rank_boundary, 2040, comm, ierr)
          WRITE(STAT%FHNDL,*) 'MPI_SEND ierr=', ierr
          FLUSH(STAT%FHNDL)
        END IF
      ELSE
        IF (rank_boundary .ne. rank_hasboundary) THEN
          IF (myrank .eq. rank_hasboundary) THEN
            CALL MPI_SEND(WBAC, MSC*MDC, rtype, rank_boundary, 2035, comm, ierr)
          END IF
          IF (myrank .eq. rank_boundary) THEN
            CALL MPI_RECV(WBAC,MSC*MDC,rtype, rank_hasboundary, 2035, comm, istatus, ierr)
          END IF
        END IF
        IF (myrank .eq. rank_boundary) THEN
          DO IP=1,IWBMNPGL
            WBAC_GL(:,:,IP)=WBAC(:,:,1)
          END DO
        END IF
      END IF
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE TRIG_SYNCHRONIZATION(V)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), intent(inout) :: V(MNEextent)
#ifdef MPI_PARALL_GRID
      integer iNeigh, iRank
      DO iNeigh=1,ie_nnbr_send
        iRank=ListNeigh_ie_send(iNeigh)
        CALL mpi_isend(V, 1, ie_send_type(iNeigh), iRank, 1020, comm, ie_send_rqst(iNeigh), ierr)
      END DO
      DO iNeigh=1,ie_nnbr_recv
        iRank=ListNeigh_ie_recv(iNeigh)
        call mpi_irecv(V,1,ie_recv_type(iNeigh),iRank,1020,comm,ie_recv_rqst(iNeigh),ierr)
      END DO
      IF (ie_nnbr_send > 0) THEN
        call mpi_waitall(ie_nnbr_send, ie_send_rqst, ie_send_stat,ierr)
      END IF
      IF (ie_nnbr_recv > 0) THEN
        call mpi_waitall(ie_nnbr_recv, ie_recv_rqst, ie_recv_stat,ierr)
      END IF
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef MPI_PARALL_GRID
# if defined NETCDF && defined DEBUG
      SUBROUTINE NETCDF_WRITE_MATRIX(LocalColor, ASPAR)
      USE DATAPOOL
      USE NETCDF
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      integer, SAVE :: iSystem = 1
      integer MSCeffect
      integer ired, ncid, var_id
      MSCeffect=LocalColor%MSCeffect
      WRITE (FILE_NAME,10) TRIM(PRE_FILE_NAME),nproc, iSystem, myrank
  10  FORMAT (a,'_np',i2.2,'_syst',i3.3,'_iproc',i4.4, '.nc')
      iret = nf90_create(TRIM(FILE_NAME), NF90_CLOBBER, ncid)
      iret = nf90_def_dim(ncid, 'iter', NF90_UNLIMITED, iter_dims)
      iret = nf90_def_dim(ncid, 'three', 3, three_dims)
      iret = nf90_def_dim(ncid, 'nfreq', MSCeffect, nfreq_dims)
      iret = nf90_def_dim(ncid, 'ndir', MDC, ndir_dims)
      iret = nf90_def_dim(ncid, 'bbz', MNP, mnp_dims)
      iret = nf90_def_dim(ncid, 'mnpp', MNP+1, mnpp_dims)
      iret = nf90_def_dim(ncid, 'np_global', np_global, npgl_dims)
      iret = nf90_def_dim(ncid, 'np_res', NP_RES, np_res_dims)
      iret = nf90_def_dim(ncid, 'mne', MNE, mne_dims)
      iret = nf90_def_dim(ncid, 'nnz', NNZ, nnz_dims)
      iret = nf90_def_var(ncid, 'ASPAR', NF90_DOUBLE,(/nfreq_dims, ndir_dims,nnz_dims/),var_id)
      iret = nf90_def_var(ncid, 'IA', NF90_INT,(/mnpp_dims/),var_id)
      iret = nf90_def_var(ncid, 'JA', NF90_INT,(/nnz_dims/),var_id)
      iret = nf90_def_var(ncid, 'iplg', NF90_INT,(/ mnp_dims/),var_id)
      iret = nf90_close(ncid)
      !
      iret = nf90_open(TRIM(FILE_NAME), NF90_WRITE, ncid)
      iret=nf90_inq_varid(ncid, 'iplg', var_id)
      iret=nf90_put_var(ncid,var_id,iplg,start=(/1/), count = (/ MNP /))
      iret=nf90_inq_varid(ncid, 'IA', var_id)
      iret=nf90_put_var(ncid,var_id,IA,start=(/1/), count = (/ MNP+1 /))
      !
      iret=nf90_inq_varid(ncid, 'JA', var_id)
      iret=nf90_put_var(ncid,var_id,JA,start=(/1/), count = (/ NNZ /))
      !
      iret=nf90_inq_varid(ncid, 'ine', var_id)
      iret=nf90_put_var(ncid,var_id,INE,start=(/1,1/), count = (/ 3, MNE /))
      !
      iret=nf90_inq_varid(ncid, 'ASPAR', var_id)
      iret=nf90_put_var(ncid,var_id,ASPAR,start=(/1,1,1/), count = (/ MSC, MDC, NNZ/))
      iret=nf90_inq_varid(ncid, 'IA', var_id)
      iret=nf90_put_var(ncid,var_id,IA,start=(/1/), count = (/ MNP+1/))
      iret=nf90_inq_varid(ncid, 'JA', var_id)
      iret=nf90_put_var(ncid,var_id,JA,start=(/1/), count = (/ NNZ/))
      iret=nf90_inq_varid(ncid, 'ListPos', var_id)
      iret=nf90_put_var(ncid,var_id,ListPos,start=(/1/), count = (/ np_global/))
      iret = nf90_close(ncid)
      !
      END SUBROUTINE
# endif
#endif
