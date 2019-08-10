#include "wwm_functions.h"
! I5 is under the assumptions that we have a lot of memory
!    and tries to minimize the number of MPI exchanges and
!    to have the processors be as busy as possible
!    We use memory ordered as AC(MNP,MSC,MDC)
! I5B is the same as I5. We use memory ordered as AC(MSC,MDC,MNP)
!    so reordering at the beginning but less operations later on.
#undef DEBUG
!#define DEBUG
! This is for the reordering of ASPAR_pc and hopefully higher speed
! in the application of the preconditioner.
#undef REORDER_ASPAR_PC
#define REORDER_ASPAR_PC
! This is for the computation of ASPAR_block by a block algorithm
! with hopefully higher speed.
#undef ASPAR_B_COMPUTE_BLOCK
#define ASPAR_B_COMPUTE_BLOCK
! Either we use the SCHISM exchange routine or ours that exchanges only
! the ghost nodes and not the interface nodes.
#undef NO_SELFE_EXCH
#define NO_SELFE_EXCH
! Repeated CX/CY computations but less memory used.
#undef NO_MEMORY_CX_CY
#define NO_MEMORY_CX_CY
! For the SOR preconditioner, we can actually compute directly
! from the matrix since it is so simple.
#undef SOR_DIRECT
#define SOR_DIRECT
!
#undef SINGLE_LOOP_AMATRIX
#define SINGLE_LOOP_AMATRIX
!
! More complexity! Some options excludes other!
!
#if defined REORDER_ASPAR_PC && defined SOR_DIRECT
# undef REORDER_ASPAR_PC
#endif
!**********************************************************************
!* We have to think on how the system is solved. Many questions are   *
!* mixed: the ordering of the nodes, the ghost nodes, the aspar array *
!* Here is a repository of the conclusions that have been reached     *
!*                                                                    *
!* Ordering 1> should be that way: We have two global nodes i and j.  *
!* -- If i and j belong to a common local grid, then we select the    *
!*    grid of lowest color and decide whether ipgl(i) < ipgl(j)       *
!* -- If i and j belong to two different grid then                    *
!*     ---If Color(i) < Color(j) or reverse we decide by that         *
!*     ---If Color(i) = Color(j) we decide by i<j or not (but it      *
!*        does not matter to the solution)                            *
!* The functions WRITE_EXPLICIT_ORDERING does exactly that and        *
!* provides an ordering that can be used. That is we start with the   *
!* nodes of lowest color and index until we arrive at highest color   *
!*                                                                    *
!* The ASPAR is computed correctly only on 1:NP_RES but this can be   *
!* extended by exchange routines.                                     *
!* We Compute on the resident nodes only. This means loops over       *
!* IP=1,NP_RES and backwards. This means that we do not have to do    *
!* exchanges of ASPAR values. Only the resident nodes are sent.       *
!* This is smaller and this is all that we ever need.                 *
!*                                                                    *
!* WRONG APPROACHES:                                                  *
!* to use all nodes 1:MNP may look simpler but it forces to have the  *
!* following property of the NNZ, IA, JA arrays. If two vertices      *
!* i and j are adjacent in a grid G, then they are adjacent in ANY    *
!* of the grid in which they are both contained.                      *
!* This property is actually not satisfied in general.                *
!* We could extend the IA, JA arrays by                               *
!* adding some vertices but that looks quite hazardous idea and it    *
!* it is actually not needed by the ILU0 preconditioner and other     *
!*                                                                    *
!* PROBLEM:                                                           *
!* There is an asymmetry in the construction of the ordering.         *
!* We start from low colors and upwards. If we had started with       *
!* high colors and gone downwards, then we get a different ordering   *
!* (even if we take the opposite, because of the interface nodes)     *
!* This requires the construction of many mappings.                   *
!* Our approach is actually to rebuild separate node sets.            *
!*                                                                    *
!* CHECKS:                                                            *
!* ---The sum of number of non-zero entries in Jstatus_L over all     *
!*    nodes should be equal to the sum of number of non-zero entries  *
!*    of Jstatus_U over all nodes.                                    *
!*    This is because number of upper diagonal entries should be      *
!*    equal to number of lower diagonal entries.                      *
!* ---We CANNOT have Jstatus_L(J)=1 and Jstatus_U(J)=1, i.e. a matrix *
!*    entry cannot be both lower and upper.                           *
!* ---We have sum of Jstatus_L + sum J_status_U + np_global should    *
!*    be equal to NNZ_global                                          *
!*                                                                    *
!* So, procedure is as follows:                                       *
!* ---compute ASPAR on 1,NP_RES nodes and no synchronization          *
!* ---compute on IP=1,NP_RES for L solving                            *
!* ---export to grids of higher rank, the values on nodes 1,NP_RES    *
!*    only. The other ghost points have invalid values or are         *
!*    resident of other grids of lower rank. (at this stage, some     *
!*    ghost values are wrong but are not exported)                    *
!* ---export to grid of lower rank in order to correct their ghost    *
!*    values and get the value                                        *
!* ---compute on IP=NP_RES,1,-1 for U solving                         *
!* ---export to grid of lower rank.                                   *
!*    Do everything similarly to L solve.                             *
!*                                                                    *
!* The basic approach is that we compute at a node S if and only if   *
!* it is a resident node. The twist come because some nodes are       *
!* resident for TWO domains. This is why we have the CovLower         *
!* We need to create disjoint domains for each node, so that          *
!* we have sum   sum(CovLower) = np_global                            *
!*                                                                    *
!* At the end of the resolution of the system, we need to do the      *
!* synchronization with respect to the unused nodes.                  *
!* Mystery?: When we apply the function, we do it on 1:NP_RES and     *
!* then call synchronizer. So, this means we need to do the           *
!* synchronization after the call to the preconditioner.              *
!* But there may be space for improvements here.                      *
!*                                                                    *
!* Description of specific exchange arrays:                           *
!* ---wwm_p2dsend_type/wwm_p2drecv_type                               *
!*    The points of 1:NP_RES are sent to nodes that contained them    *
!*    length=1                                                        *
!* ---wwmtot_p2dsend_type/wwmtot_p2drecv_type                         *
!*    same as above but length=MSC*MDC                                *
!* ---blk_p2dsend_type/blk_p2drecv_type                               *
!*    same as above but length=maxBlockLength for matrix exchanges    *
!* ---wwmmat_p2dsend_type/wwmmat_p2drecv_type                         *
!*    exchange of correct matrix elements, i.e. elements A(I,J)       *
!*    with I<=NP_RES                                                  *
!*    length is 1.                                                    *
!* ---u2l_p2dsend_type/u2l_p2drecv_type                               *
!*    upper 2 lower exchange arrays, depends on CovLower, so on       *
!*    the coloring chosen. length=maxBlockLength                      *
!*    exchange are from upper to lower.                               *
!* ---sync_p2dsend_type/sync_p2drecv_type                             *
!*    synchronize value, i.e. the CovLower=1 values  are send to all  *
!*    nodes. Length is MSC*MDC                                        *
!**********************************************************************
#if defined WWM_SOLVER && defined MPI_PARALL_GRID
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_EXCHANGE_P4D_WWM(LocalColor, AC)
      USE DATAPOOL, only : MSC, MDC, rkind, LocalColorInfo
      USE DATAPOOL, only : wwm_nnbr_send, wwm_nnbr_recv
      USE DATAPOOL, only : wwm_ListNeigh_send, wwm_ListNeigh_recv
      USE DATAPOOL, only : wwmtot_p2dsend_type, wwmtot_p2drecv_type
      USE DATAPOOL, only : wwm_p2dsend_type, wwm_p2drecv_type
      USE DATAPOOL, only : wwm_p2dsend_rqst, wwm_p2drecv_rqst
      USE DATAPOOL, only : wwm_p2dsend_stat, wwm_p2drecv_stat
      USE DATAPOOL, only : ZERO, NP_RES, MNP
      USE datapool, only : comm, ierr, myrank
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      real(rkind), intent(inout) :: AC(MSC,MDC,MNP)
      integer iSync, iRank
      integer IS, ID, IP
      real(rkind) SumErr, Lerror
      real(rkind) :: ACtest(MSC,MDC,MNP)
      real(rkind) :: U1(MNP), U2(MNP)
      CALL I5B_TOTAL_COHERENCY_ERROR_NPRES(MSC, AC, Lerror)
!      Print *, 'NP_RES cohenrency error=', Lerror
#ifdef DEBUG
      WRITE(740+myrank,*) 'I5B_EXCHANGE_P4D_WWM, begin, sum(AC)=', sum(AC)
      FLUSH(740+myrank)
#endif
      ACtest=AC
      SumErr=ZERO
      DO IS=1,MSC
        DO ID=1,MDC
          U1=AC(IS,ID,:)
          U2=AC(IS,ID,:)
          DO iSync=1,wwm_nnbr_send
            iRank=wwm_ListNeigh_send(iSync)
            CALL mpi_isend(U1, 1, wwm_p2dsend_type(iSync), iRank-1, 1020, comm, wwm_p2dsend_rqst(iSync), ierr)
          END DO
          DO iSync=1,wwm_nnbr_recv
            iRank=wwm_ListNeigh_recv(iSync)
            call mpi_irecv(U2,1,wwm_p2drecv_type(iSync),iRank-1,1020,comm,wwm_p2drecv_rqst(iSync),ierr)
          END DO
          IF (wwm_nnbr_send > 0) THEN
            call mpi_waitall(wwm_nnbr_send, wwm_p2dsend_rqst, wwm_p2dsend_stat,ierr)
          END IF
          IF (wwm_nnbr_recv > 0) THEN
            call mpi_waitall(wwm_nnbr_recv, wwm_p2drecv_rqst, wwm_p2drecv_stat,ierr)
          END IF
          DO IP=1,NP_RES
            SumErr=SumErr + abs(U1(IP) - U2(IP))
          END DO
        END DO
      END DO
#ifdef DEBUG
      WRITE(740+myrank,*) 'Total SumErr=', SumErr
      FLUSH(740+myrank)
#endif
      DO iSync=1,wwm_nnbr_send
        iRank=wwm_ListNeigh_send(iSync)
        CALL mpi_isend(ACtest, 1, wwmtot_p2dsend_type(iSync), iRank-1, 1020, comm, wwm_p2dsend_rqst(iSync), ierr)
      END DO
      DO iSync=1,wwm_nnbr_recv
        iRank=wwm_ListNeigh_recv(iSync)
        call mpi_irecv(AC,1,wwmtot_p2drecv_type(iSync),iRank-1,1020,comm,wwm_p2drecv_rqst(iSync),ierr)
      END DO
      IF (wwm_nnbr_send > 0) THEN
        call mpi_waitall(wwm_nnbr_send, wwm_p2dsend_rqst, wwm_p2dsend_stat,ierr)
      END IF
      IF (wwm_nnbr_recv > 0) THEN
        call mpi_waitall(wwm_nnbr_recv, wwm_p2drecv_rqst, wwm_p2drecv_stat,ierr)
      END IF
      SumErr=ZERO
      DO IS=1,MSC
        DO ID=1,MDC
          DO IP=1,NP_RES
            SumErr=SumErr + abs(AC(IS,ID,IP) - ACtest(IS,ID,IP))
          END DO
        END DO
      END DO
#ifdef DEBUG
      WRITE(740+myrank,*) 'SumErr=', SumErr
      FLUSH(740+myrank)
#endif
#ifdef DEBUG
      WRITE(740+myrank,*) 'I5B_EXCHANGE_P4D_WWM, end, sum(AC)=', sum(AC)
      FLUSH(740+myrank)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_EXCHANGE_SL_WWM(LocalColor, AC)
      USE DATAPOOL, only : MSC, MDC, rkind, LocalColorInfo
      USE DATAPOOL, only : wwm_nnbr_send_sl, wwm_nnbr_recv_sl
      USE DATAPOOL, only : wwm_ListNeigh_send_sl, wwm_ListNeigh_recv_sl, wwm_ListNeigh_send
      USE DATAPOOL, only : wwmsl_send_type, wwmsl_recv_type
      USE DATAPOOL, only : wwmsl_send_rqst, wwmsl_recv_rqst
      USE DATAPOOL, only : wwmsl_send_stat, wwmsl_recv_stat
      USE DATAPOOL, only : ZERO, NP_RES, MNP
      USE datapool, only : comm, ierr, myrank
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      real(rkind), intent(inout) :: AC(MSC,MDC,MNP)
      integer iSync, iRank
      DO iSync=1,wwm_nnbr_send_sl
        iRank=wwm_ListNeigh_send(iSync)
        CALL mpi_isend(AC, 1, wwmsl_send_type(iSync), iRank-1, 1020, comm, wwmsl_send_rqst(iSync), ierr)
      END DO
      DO iSync=1,wwm_nnbr_recv_sl
        iRank=wwm_ListNeigh_recv_sl(iSync)
        call mpi_irecv(AC,1,wwmsl_recv_type(iSync),iRank-1,1020,comm,wwmsl_recv_rqst(iSync),ierr)
      END DO
      IF (wwm_nnbr_send_sl > 0) THEN
        call mpi_waitall(wwm_nnbr_send_sl, wwmsl_send_rqst, wwmsl_send_stat,ierr)
      END IF
      IF (wwm_nnbr_recv_sl > 0) THEN
        call mpi_waitall(wwm_nnbr_recv_sl, wwmsl_recv_rqst, wwmsl_recv_stat,ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_EXCHANGE_ASPAR(LocalColor, ASPAR_bl)
      USE DATAPOOL, only: comm, ierr, myrank, ierr, wwm_nnbr_m_send, wwm_ListNeigh_m_send
      use datapool, only: wwmmat_p2dsend_type, wwmmat_p2dsend_rqst, wwm_nnbr_m_recv, wwm_ListNeigh_m_recv
      use datapool, only: wwmmat_p2drecv_type, wwmmat_p2drecv_rqst, LocalColorInfo, rkind
      use datapool, only: wwmmat_p2drecv_stat, wwmmat_p2dsend_stat, MSC, MDC, NNZ
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      real(rkind), intent(inout) :: ASPAR_bl(MSC,MDC,NNZ)
      integer I, iProc
      do I=1,wwm_nnbr_m_send
        iProc=wwm_ListNeigh_m_send(I)
        call mpi_isend(ASPAR_bl,1,wwmmat_p2dsend_type(I),iProc,991,comm,wwmmat_p2dsend_rqst(i),ierr)
      enddo
      do I=1,wwm_nnbr_m_recv
        iProc=wwm_ListNeigh_m_recv(I)
        call mpi_irecv(ASPAR_bl,1,wwmmat_p2drecv_type(I),iProc,991,comm,wwmmat_p2drecv_rqst(i),ierr)
      enddo
      IF (wwm_nnbr_m_recv .gt. 0) THEN
        call mpi_waitall(wwm_nnbr_m_recv,wwmmat_p2drecv_rqst,wwmmat_p2drecv_stat,ierr)
      END IF
      IF (wwm_nnbr_m_send .gt. 0) THEN
        call mpi_waitall(wwm_nnbr_m_send,wwmmat_p2dsend_rqst,wwmmat_p2dsend_stat,ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COLLECT_ALL_COVLOWER(LocalColor)
      USE DATAPOOL
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      integer, allocatable :: rbuf_int(:)
      integer len, iProc, IP, idx, sumMNP
      sumMNP=sum(ListMNP)
      allocate(LocalColor % ListCovLower(sumMNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 18')
      IF (myrank == 0) THEN
        idx=0
        DO IP=1,MNP
          idx=idx+1
          LocalColor % ListCovLower(idx)=LocalColor % CovLower(IP)
        END DO
        DO iProc=2,nproc
          len=ListMNP(iProc)
          allocate(rbuf_int(len), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 19')
          CALL MPI_RECV(rbuf_int,len,itype, iProc-1, 809, comm, istatus, ierr)
          DO IP=1,len
            idx=idx+1
            LocalColor % ListCovLower(idx)=rbuf_int(IP)
          END DO
          deallocate(rbuf_int)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(LocalColor % ListCovLower,sumMNP,itype, iProc-1, 811, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(LocalColor % CovLower,MNP,itype, 0, 809, comm, ierr)
        CALL MPI_RECV(LocalColor % ListCovLower,sumMNP,itype, 0, 811, comm, istatus, ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CREATE_WWM_P2D_EXCH
      USE DATAPOOL
      implicit none
      integer :: ListFirst(nproc)
      integer :: ListCommon_recv(nproc)
      integer :: ListCommon_send(nproc)
      integer :: ListMapped(np_global)
      integer :: ListMappedB(np_global)
      integer, allocatable :: dspl_send(:), dspl_recv(:)
      integer, allocatable :: dspl_send_tot(:), dspl_recv_tot(:)
      integer IP, IP_glob, iProc, MNPloc, idx, NP_RESloc
      integer iNeigh, IPmap, nbCommon
      integer nbCommon_send, nbCommon_recv, idx_send, idx_recv
      integer nbCommon_send_sl, nbCommon_recv_sl
      integer sumNbCommon_send, sumNbCommon_recv
      integer idxDspl_send, idxDspl_recv
      ListFirst=0
      DO iProc=2,nproc
        ListFirst(iProc)=ListFirst(iProc-1) + ListMNP(iProc-1)
      END DO
      !
      ! First the arrays so that coordinates 1:NP_RES got send 
      ! to all nodes 1:MNP of neighboring components
      !
      ListMapped=0
      DO IP=1,MNP
        IP_glob=iplg(IP)
        ListMapped(IP_glob)=IP
      END DO
      wwm_nnbr_send=0
      wwm_nnbr_recv=0
      ListCommon_send=0
      ListCommon_recv=0
      DO iProc=1,nproc
        IF (iPROC .ne. myrank+1) THEN
          MNPloc=ListMNP(iProc)
          NP_RESloc=ListNP_RES(iProc)
          ListMappedB=0
          DO IP=1,MNPloc
            IP_glob=ListIPLG(IP+ListFirst(iProc))
            ListMappedB(IP_glob)=IP
          END DO
          !
          nbCommon_recv=0
          DO IP=1,NP_RESloc
            IP_glob=ListIPLG(IP+ListFirst(iProc))
            IF (ListMapped(IP_glob).gt.0) THEN
              nbCommon_recv=nbCommon_recv+1
            END IF
          END DO
          IF (nbCommon_recv .gt. 0) THEN
            wwm_nnbr_recv=wwm_nnbr_recv+1
            ListCommon_recv(iProc)=nbCommon_recv
          END IF
          !
          nbCommon_send=0
          DO IP=1,NP_RES
            IP_glob=iplg(IP)
            IF (ListMappedB(IP_glob).gt.0) THEN
              nbCommon_send=nbCommon_send+1
            END IF
          END DO
          IF (nbCommon_send .gt. 0) THEN
            wwm_nnbr_send=wwm_nnbr_send+1
            ListCommon_send(iProc)=nbCommon_send
          END IF
        END IF
      END DO
      allocate(wwm_ListNbCommon_send(wwm_nnbr_send), wwm_ListNbCommon_recv(wwm_nnbr_recv), wwm_ListNeigh_send(wwm_nnbr_send), wwm_ListNeigh_recv(wwm_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 27')
      idx_send=0
      idx_recv=0
      sumNbCommon_send=0
      sumNbCommon_recv=0
      DO iProc=1,nproc
        IF (ListCommon_send(iProc) .gt. 0) THEN
          idx_send=idx_send+1
          wwm_ListNeigh_send(idx_send)=iProc
          nbCommon=ListCommon_send(iProc)
          wwm_ListNbCommon_send(idx_send)=nbCommon
          sumNbCommon_send=sumNbCommon_send+nbCommon
        END IF
        IF (ListCommon_recv(iProc) .gt. 0) THEN
          idx_recv=idx_recv+1
          wwm_ListNeigh_recv(idx_recv)=iProc
          nbCommon=ListCommon_recv(iProc)
          wwm_ListNbCommon_recv(idx_recv)=nbCommon
          sumNbCommon_recv=sumNbCommon_recv+nbCommon
        END IF
      END DO
      allocate(wwm_ListDspl_send(sumNbCommon_send), wwm_ListDspl_recv(sumNbCommon_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 28')
      wwm_nnbr=0
      DO iProc=1,nproc
        IF ((ListCommon_send(iProc).gt.0).or.(ListCommon_recv(iProc).gt.0)) THEN
          wwm_nnbr=wwm_nnbr+1
        END IF
      END DO
      allocate(wwm_ListNeigh(wwm_nnbr), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 29')
      idx=0
      DO iProc=1,nproc
        IF ((ListCommon_send(iProc).gt.0).or.(ListCommon_recv(iProc).gt.0)) THEN
          idx=idx+1
          wwm_ListNeigh(idx)=iProc
        END IF
      END DO
      allocate(wwm_p2dsend_rqst(wwm_nnbr_send), wwm_p2drecv_rqst(wwm_nnbr_recv), &
     &wwm_p2dsend_stat(MPI_STATUS_SIZE,wwm_nnbr_send), wwm_p2drecv_stat(MPI_STATUS_SIZE,wwm_nnbr_recv), &
     &wwm_p2dsend_type(wwm_nnbr_send), wwm_p2drecv_type(wwm_nnbr_recv), &
     &wwmtot_p2dsend_type(wwm_nnbr_send), wwmtot_p2drecv_type(wwm_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 30')
      idxDspl_send=0
      DO iNeigh=1,wwm_nnbr_send
        iProc=wwm_ListNeigh_send(iNeigh)
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        ListMappedB=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          ListMappedB(IP_glob)=IP
        END DO
        nbCommon=wwm_ListNbCommon_send(iNeigh)
        allocate(dspl_send(nbCommon), dspl_send_tot(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 31')
        idx=0
        DO IP=1,NP_RES
          IP_glob=iplg(IP)
          IF (ListMappedB(IP_glob).gt.0) THEN
            IPmap=ListMappedB(IP_glob)
            idx=idx+1
            dspl_send(idx)=IP-1
            dspl_send_tot(idx)=MSC*MDC*(IP-1)
            idxDspl_send=idxDspl_send+1
            wwm_ListDspl_send(idxDspl_send)=IP
          END IF
        END DO
        call mpi_type_create_indexed_block(nbCommon,1,dspl_send,rtype,wwm_p2dsend_type(iNeigh), ierr)
        call mpi_type_commit(wwm_p2dsend_type(iNeigh), ierr)
        call mpi_type_create_indexed_block(nbCommon,MSC*MDC,dspl_send_tot,rtype,wwmtot_p2dsend_type(iNeigh), ierr)
        call mpi_type_commit(wwmtot_p2dsend_type(iNeigh), ierr)
        deallocate(dspl_send, dspl_send_tot)
      END DO
      idxDspl_recv=0
      DO iNeigh=1,wwm_nnbr_recv
        iProc=wwm_ListNeigh_recv(iNeigh)
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        nbCommon=wwm_ListNbCommon_recv(iNeigh)
        allocate(dspl_recv(nbCommon), dspl_recv_tot(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 32')
        idx=0
        DO IP=1,NP_RESloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          IF (ListMapped(IP_glob).gt.0) THEN
            IPmap=ListMapped(IP_glob)
            idx=idx+1
            dspl_recv(idx)=IPmap-1
            dspl_recv_tot(idx)=MSC*MDC*(IPmap-1)
            idxDspl_recv=idxDspl_recv+1
            wwm_ListDspl_recv(idxDspl_recv)=IPmap
#ifdef DEBUG
            WRITE(800+myrank,*) 'Recv IP=', IPmap, 'IPglob=', IP_glob
            FLUSH(800+myrank)
#endif 
          END IF
        END DO
        !
        call mpi_type_create_indexed_block(nbCommon,1,dspl_recv,rtype,wwm_p2drecv_type(iNeigh),ierr)
        call mpi_type_commit(wwm_p2drecv_type(iNeigh), ierr)
        call mpi_type_create_indexed_block(nbCommon,MSC*MDC,dspl_recv_tot,rtype,wwmtot_p2drecv_type(iNeigh),ierr)
        call mpi_type_commit(wwmtot_p2drecv_type(iNeigh), ierr)
        deallocate(dspl_recv, dspl_recv_tot)
      END DO
      !
      ! First the arrays so that coordinates 1:NP_RES got send 
      ! to all nodes NP_RES+1:MNP of neighboring components
      ! "sl" : "super local"
      !
      ListMapped=0
      DO IP=NP_RES+1,MNP
        IP_glob=iplg(IP)
        ListMapped(IP_glob)=IP
      END DO
      wwm_nnbr_send_sl=0
      wwm_nnbr_recv_sl=0
      ListCommon_send=0
      ListCommon_recv=0
      DO iNeigh=1,wwm_nnbr
        iProc=wwm_ListNeigh(iNeigh)
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        ListMappedB=0
        DO IP=NP_RESloc+1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          ListMappedB(IP_glob)=IP
        END DO
        !
        nbCommon_recv_sl=0
        DO IP=1,NP_RESloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          IF (ListMapped(IP_glob).gt.0) THEN
            nbCommon_recv_sl=nbCommon_recv_sl+1
          END IF
        END DO
        IF (nbCommon_recv_sl .gt. 0) THEN
          wwm_nnbr_recv_sl=wwm_nnbr_recv_sl+1
          ListCommon_recv(iProc)=nbCommon_recv_sl
        END IF
        !
        nbCommon_send_sl=0
        DO IP=1,NP_RES
          IP_glob=iplg(IP)
          IF (ListMappedB(IP_glob).gt.0) THEN
            nbCommon_send_sl=nbCommon_send_sl+1
          END IF
        END DO
        IF (nbCommon_send_sl .gt. 0) THEN
          wwm_nnbr_send_sl=wwm_nnbr_send_sl+1
          ListCommon_send(iProc)=nbCommon_send_sl
        END IF
      END DO
      allocate(wwm_ListNbCommon_send_sl(wwm_nnbr_send_sl), wwm_ListNbCommon_recv_sl(wwm_nnbr_recv_sl), wwm_ListNeigh_send_sl(wwm_nnbr_send_sl), wwm_ListNeigh_recv_sl(wwm_nnbr_recv_sl), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 33')
      idx_send=0
      idx_recv=0
      DO iProc=1,nproc
        nbCommon=ListCommon_send(iProc)
        IF (nbCommon .gt. 0) THEN
          idx_send=idx_send+1
          wwm_ListNeigh_send_sl(idx_send)=iProc
          wwm_ListNbCommon_send_sl(idx_send)=nbCommon
        END IF
        nbCommon=ListCommon_recv(iProc)
        IF (nbCommon .gt. 0) THEN
          idx_recv=idx_recv+1
          wwm_ListNeigh_recv_sl(idx_recv)=iProc
          wwm_ListNbCommon_recv_sl(idx_recv)=nbCommon
        END IF
      END DO
      allocate(wwmsl_send_rqst(wwm_nnbr_send_sl), wwmsl_recv_rqst(wwm_nnbr_recv_sl), &
     &wwmsl_send_stat(MPI_STATUS_SIZE,wwm_nnbr_send_sl), wwmsl_recv_stat(MPI_STATUS_SIZE,wwm_nnbr_recv_sl), &
     &wwmsl_send_type(wwm_nnbr_send_sl), wwmsl_recv_type(wwm_nnbr_recv_sl), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 35')
      idxDspl_send=0
      DO iNeigh=1,wwm_nnbr_send_sl
        iProc=wwm_ListNeigh_send_sl(iNeigh)
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        ListMappedB=0
        DO IP=NP_RESloc+1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          ListMappedB(IP_glob)=IP
        END DO
        nbCommon=wwm_ListNbCommon_send_sl(iNeigh)
        allocate(dspl_send(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 36')
        idx=0
        DO IP=1,NP_RES
          IP_glob=iplg(IP)
          IF (ListMappedB(IP_glob).gt.0) THEN
            IPmap=ListMappedB(IP_glob)
            idx=idx+1
            dspl_send(idx)=MSC*MDC*(IP-1)
          END IF
        END DO
        call mpi_type_create_indexed_block(nbCommon,MSC*MDC,dspl_send,rtype,wwmsl_send_type(iNeigh), ierr)
        call mpi_type_commit(wwmsl_send_type(iNeigh), ierr)
        deallocate(dspl_send)
      END DO
      idxDspl_recv=0
      DO iNeigh=1,wwm_nnbr_recv_sl
        iProc=wwm_ListNeigh_recv_sl(iNeigh)
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        nbCommon=wwm_ListNbCommon_recv_sl(iNeigh)
        allocate(dspl_recv(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 37')
        idx=0
        DO IP=1,NP_RESloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          IF (ListMapped(IP_glob).gt.0) THEN
            IPmap=ListMapped(IP_glob)
            idx=idx+1
            dspl_recv(idx)=MSC*MDC*(IPmap-1)
          END IF
        END DO
        call mpi_type_create_indexed_block(nbCommon,MSC*MDC,dspl_recv,rtype,wwmsl_recv_type(iNeigh),ierr)
        call mpi_type_commit(wwmsl_recv_type(iNeigh), ierr)
        deallocate(dspl_recv)
      END DO
      END SUBROUTINE
!**********************************************************************
!*  This subroutine creates array for exchanging matrix entries       *
!*  Process cannot get all the matrix entries right, and so they need *
!*  to exchange data.                                                 *
!*  Matrix entries A(i,j) with i <= NP_RES are correct and are        *
!*  exported                                                          *
!**********************************************************************
      SUBROUTINE CREATE_WWM_MAT_P2D_EXCH
      USE DATAPOOL
      implicit none
      integer :: ListFirstMNP(nproc), ListFirstNNZ(nproc)
      integer :: ListCommon_send(nproc), ListCommon_recv(nproc)
      integer :: ListMapped(np_global)
      integer :: ListMappedB(np_global)
      integer, allocatable :: dspl_send(:), dspl_recv(:)
      integer nbCommon_send, nbCommon_recv
      integer IAfirst
      integer IP, JP, I, J, J2, IP_glob, JP_glob, iProc
      integer MNPloc, NP_RESloc, JP_j
      integer IPloc, JPloc, Jfound, idx
      integer iNeigh, nbCommon
      integer sumNbCommon_send, sumNbCommon_recv
      integer idxDspl_send, idxDspl_recv
      ListFirstNNZ=0
      ListFirstMNP=0
      DO iProc=2,nproc
        ListFirstMNP(iProc)=ListFirstMNP(iProc-1) + ListMNP(iProc-1)
        ListFirstNNZ(iProc)=ListFirstNNZ(iProc-1) + ListNNZ(iProc-1)
      END DO
      ListMapped=0
      DO IP=1,MNP
        IP_glob=iplg(IP)
        ListMapped(IP_glob)=IP
      END DO
      wwm_nnbr_m_recv=0
      wwm_nnbr_m_send=0
      ListCommon_recv=0
      DO I=1,wwm_nnbr_recv
        iProc=wwm_ListNeigh_recv(I)
        NP_RESloc=ListNP_RES(iProc)
        MNPloc=ListMNP(iProc)
        ListMappedB=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirstMNP(iProc))
          ListMappedB(IP_glob)=IP
        END DO
        nbCommon_recv=0
        DO IP=1,NP_RESloc
          IP_glob=ListIPLG(IP+ListFirstMNP(iProc))
          IPloc=ListMapped(IP_glob)
          IF (IPloc.gt.0) THEN
            IAfirst=ListFirstMNP(iProc) + iProc-1
            DO J=ListIA(IP+IAfirst),ListIA(IP+IAfirst+1)-1
              JP=ListJA(J+ListFirstNNZ(iProc))
              JP_glob=ListIPLG(JP+ListFirstMNP(iProc))
              JPloc=ListMapped(JP_glob)
              IF (JPloc.gt.0) THEN
                JFOUND=-1
                DO J2=IA(IPloc),IA(IPloc+1)-1
                  IF (JA(J2) == JPloc) THEN
                    JFOUND=J2
                  END IF
                END DO
                IF (JFOUND .gt. 0) THEN
                  nbCommon_recv=nbCommon_recv+1
                END IF
              END IF
            END DO
          END IF
        END DO
        IF (nbCommon_recv .gt. 0) THEN
          wwm_nnbr_m_recv=wwm_nnbr_m_recv+1
          ListCommon_recv(iProc)=nbCommon_recv
        END IF
      END DO
      ListCommon_send=0
      DO I=1,wwm_nnbr_send
        iProc=wwm_ListNeigh_send(I)
        NP_RESloc=ListNP_RES(iProc)
        MNPloc=ListMNP(iProc)
        ListMappedB=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirstMNP(iProc))
          ListMappedB(IP_glob)=IP
        END DO
        nbCommon_send=0
        DO IP=1,NP_RES
          IP_glob=iplg(IP)
          IPloc=ListMappedB(IP_glob)
          IF (IPloc.gt.0) THEN
            DO J=IA(IP),IA(IP+1)-1
              JP=JA(J)
              JP_glob=iplg(JP)
              JPloc=ListMappedB(JP_glob)
              IF (JPloc.gt.0) THEN
                IAfirst=ListFirstMNP(iProc) + iProc-1
                JFOUND=-1
                DO J2=ListIA(IPloc+IAfirst),ListIA(IPloc+IAfirst+1)-1
                  JP_j=ListJA(J2+ListFirstNNZ(iProc))
                  IF (JP_j == JPloc) THEN
                    JFOUND=J2
                  END IF
                END DO
                IF (JFOUND .gt. 0) THEN
                  nbCommon_send=nbCommon_send+1
                END IF
              END IF
            END DO
          END IF
        END DO
        IF (nbCommon_send .gt. 0) THEN
          wwm_nnbr_m_send=wwm_nnbr_m_send+1
          ListCommon_send(iProc)=nbCommon_send
        END IF
      END DO
      allocate(wwmmat_p2dsend_rqst(wwm_nnbr_m_send), wwmmat_p2drecv_rqst(wwm_nnbr_m_recv), &
     &wwmmat_p2dsend_stat(MPI_STATUS_SIZE,wwm_nnbr_m_send), &
     &wwmmat_p2drecv_stat(MPI_STATUS_SIZE,wwm_nnbr_m_recv), wwmmat_p2dsend_type(wwm_nnbr_m_send), &
     &wwmmat_p2drecv_type(wwm_nnbr_m_recv), wwm_ListNbCommon_m_send(wwm_nnbr_m_send), &
     &wwm_ListNbCommon_m_recv(wwm_nnbr_m_recv), wwm_ListNeigh_m_recv(wwm_nnbr_m_recv), &
     &wwm_ListNeigh_m_send(wwm_nnbr_m_send), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 38')
      idx=0
      sumNbCommon_send=0
      DO iProc=1,nproc
        IF (ListCommon_send(iProc) .gt. 0) THEN
          idx=idx+1
          wwm_ListNeigh_m_send(idx)=iProc-1
          nbCommon=ListCommon_send(iProc)
          wwm_ListNbCommon_m_send(idx)=nbCommon
          sumNbCommon_send=sumNbCommon_send+nbCommon
        END IF
      END DO
      idx=0
      sumNbCommon_recv=0
      DO iProc=1,nproc
        IF (ListCommon_recv(iProc) .gt. 0) THEN
          idx=idx+1
          wwm_ListNeigh_m_recv(idx)=iProc-1
          nbCommon=ListCommon_recv(iProc)
          wwm_ListNbCommon_m_recv(idx)=nbCommon
          sumNbCommon_recv=sumNbCommon_recv+nbCommon
        END IF
      END DO
      allocate(wwm_ListDspl_m_send(sumNbCommon_send), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 39')
      idxDspl_send=0
      DO iNeigh=1,wwm_nnbr_m_send
        iProc=wwm_ListNeigh_m_send(iNeigh)+1
        MNPloc=ListMNP(iProc)
        ListMappedB=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirstMNP(iProc))
          ListMappedB(IP_glob)=IP
        END DO
        nbCommon=wwm_ListNbCommon_m_send(iNeigh)
        allocate(dspl_send(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 40')
        idx=0
        DO IP=1,NP_RES
          IP_glob=iplg(IP)
          IPloc=ListMappedB(IP_glob)
          IF (IPloc.gt.0) THEN
            DO J=IA(IP),IA(IP+1)-1
              JP=JA(J)
              JP_glob=iplg(JP)
              JPloc=ListMappedB(JP_glob)
              IF (JPloc .gt. 0) THEN
                IAfirst=ListFirstMNP(iProc) + iProc-1
                JFOUND=-1
                DO J2=ListIA(IPloc+IAfirst),ListIA(IPloc+IAfirst+1)-1
                  JP_j=ListJA(J2+ListFirstNNZ(iProc))
                  IF (JP_j == JPloc) THEN
                    JFOUND=J2
                  END IF
                END DO
                IF (JFOUND .gt. 0) THEN
                  idxDspl_send=idxDspl_send+1
                  wwm_ListDspl_m_send(idxDspl_send)=J
                  idx=idx+1
                  dspl_send(idx)=(J-1)*MSC*MDC
                END IF
              END IF
            END DO
          END IF
        END DO
        call mpi_type_create_indexed_block(nbCommon,MSC*MDC,dspl_send,rtype,wwmmat_p2dsend_type(iNeigh), ierr)
        call mpi_type_commit(wwmmat_p2dsend_type(iNeigh), ierr)
        deallocate(dspl_send)
      END DO
      allocate(wwm_ListDspl_m_recv(sumNbCommon_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 41')
      idxDspl_recv=0
      DO iNeigh=1,wwm_nnbr_m_recv
        iProc=wwm_ListNeigh_m_recv(iNeigh)+1
        nbCommon=wwm_ListNbCommon_m_recv(iNeigh)
        allocate(dspl_recv(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 42')
        NP_RESloc=ListNP_RES(iProc)
        idx=0
        DO IP=1,NP_RESloc
          IP_glob=ListIPLG(IP+ListFirstMNP(iProc))
          IPloc=ListMapped(IP_glob)
          IF (IPloc.gt.0) THEN
            IAfirst=ListFirstMNP(iProc) + iProc-1
            DO J=ListIA(IP+IAfirst),ListIA(IP+IAfirst+1)-1
              JP=ListJA(J+ListFirstNNZ(iProc))
              JP_glob=ListIPLG(JP+ListFirstMNP(iProc))
              JPloc=ListMapped(JP_glob)
              IF (JPloc.gt.0) THEN
                JFOUND=-1
                DO J2=IA(IPloc),IA(IPloc+1)-1
                  IF (JA(J2) == JPloc) THEN
                    JFOUND=J2
                  END IF
                END DO
                IF (JFOUND /= -1) THEN
                  idxDspl_recv=idxDspl_recv+1
                  wwm_ListDspl_m_recv(idxDspl_recv)=Jfound
                  idx=idx+1
                  dspl_recv(idx)=(Jfound-1)*MSC*MDC
                END IF
              END IF
            END DO
          END IF
        END DO
        call mpi_type_create_indexed_block(nbCommon,MSC*MDC,dspl_recv,rtype,wwmmat_p2drecv_type(iNeigh), ierr)
        call mpi_type_commit(wwmmat_p2drecv_type(iNeigh), ierr)
        deallocate(dspl_recv)
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE BUILD_MULTICOLORING(AdjGraph, ListColor)
      USE datapool, only : myrank, Graph
      implicit none
      type(Graph), intent(in) :: AdjGraph
      integer, intent(out) :: ListColor(AdjGraph%nbVert)
      integer, allocatable :: CurrColor(:)
      integer MaxDeg, iVert, eVert, eColor, eDeg
      integer idx, I, ChromaticNr, nbVert
      integer, allocatable :: ListPosFirst(:)
      integer eColorF, iVertFound, eAdjColor
      integer nbUndef, MinDeg, eAdj, MinUndef, PosMin
      integer istat
      MaxDeg=AdjGraph % MaxDeg
      nbVert=AdjGraph % nbVert
      allocate(CurrColor(MaxDeg+1), ListPosFirst(nbVert), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 43')
      ListColor=0
      idx=0
      DO iVert=1,nbVert
        ListPosFirst(iVert)=idx
        idx=idx+AdjGraph % ListDegree(iVert)
      END DO
      MinDeg=MaxDeg+3
      PosMin=-1
      DO iVert=1,nbVert
        eDeg=AdjGraph % ListDegree(iVert)
        IF (eDeg .lt. MinDeg) THEN
          MinDeg=eDeg
          PosMin=iVert
        END IF
      END DO
      idx=ListPosFirst(PosMin)
      DO I=0,MinDeg
        IF (I.eq.0) THEN
          eVert=PosMin
        ELSE
          eVert=AdjGraph % ListEdge(idx+I,2)
        END IF
        ListColor(eVert)=I+1
      END DO
      DO
        MinUndef=nbVert
        iVertFound=0
        DO iVert=1,nbVert
          IF (ListColor(iVert) == 0) THEN
            idx=ListPosFirst(iVert)
            eDeg=AdjGraph % ListDegree(iVert)
            nbUndef=0
            DO I=1,eDeg
              eAdj=AdjGraph % ListEdge(idx+I,2)
              eAdjColor=ListColor(eAdj)
              IF (eAdjColor == 0) THEN
                nbUndef=nbUndef+1
              END IF
            END DO
            IF (nbUndef .lt. MinUndef) THEN
              MinUndef=nbUndef
              iVertFound=iVert
            END IF
          END IF
        END DO
        IF (iVertFound == 0) THEN
          EXIT
        END IF
        eDeg=AdjGraph % ListDegree(iVertFound)
        idx=ListPosFirst(iVertFound)
        CurrColor=0
        DO I=1,eDeg
          eVert=AdjGraph % ListEdge(idx+I,2)
          eColor=ListColor(eVert)
          IF (eColor.gt.0) THEN
            CurrColor(eColor)=1
          END IF
        END DO
        eColorF=-1
        DO I=1,MaxDeg+1
          IF (eColorF == -1) THEN
            IF (CurrColor(I) == 0) THEN
              eColorF=I
            END IF
          END IF
        END DO
        ListColor(iVertFound)=eColorF
      END DO
      deallocate(ListPosFirst, CurrColor)
      ChromaticNr=maxval(ListColor)
# ifdef DEBUG
      WRITE(740+myrank,*) 'ChromaticNr=', ChromaticNr
      DO iVert=1,nbVert
        WRITE(740+myrank,*) 'iVert=', iVert, 'eColor=', ListColor(iVert)
      END DO
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DeallocateGraph(TheGraph)
      USE DATAPOOL, only : Graph
      implicit none
      type(Graph), intent(inout) :: TheGraph
      deallocate(TheGraph % ListDegree)
      deallocate(TheGraph % ListEdge)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_BLOCK_FREQDIR(LocalColor, Nblock)
      USE DATAPOOL, only : MNP, MDC, MSC, LocalColorInfo, rkind, stat
      USE DATAPOOL, only : wwm_nnbr_send, wwm_nnbr_recv
      USE DATAPOOL, only : wwm_ListNbCommon_send, wwm_ListNbCommon_recv
      USE DATAPOOL, only : XP, YP, rtype, ierr, myrank, iplg
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      integer, intent(in) :: Nblock
      integer Ntot, Hlen, Delta, iBlock, idx, ID, IS
      integer lenBlock, maxBlockLength
      integer istat

      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'ENTERING INIT_BLOCK_FREQDIR'
      FLUSH(STAT%FHNDL)

      Ntot=MSC*MDC
      Hlen=INT(MyREAL(Ntot)/Nblock)
      Delta=Ntot - Hlen*Nblock
      iBlock=1
      idx=1
      LocalColor % Nblock=Nblock
      IF (Delta == 0) THEN
        maxBlockLength=Hlen
      ELSE
        maxBlockLength=Hlen+1
      ENDIF
      allocate(LocalColor % ISindex(Nblock, maxBlockLength), LocalColor % IDindex(Nblock, maxBlockLength), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 44')
      DO IS=1,MSC
        DO ID=1,MDC
          LocalColor % ISindex(iBlock, idx)=IS
          LocalColor % IDindex(iBlock, idx)=ID
          IF (iBlock <= Delta) THEN
            lenBlock=Hlen+1
          ELSE
            lenBlock=Hlen
          END IF
          idx=idx+1
          IF (idx > lenBlock) THEN
            iBlock=iBlock+1
            idx=1
          ENDIF
        END DO
      END DO
      allocate(LocalColor % BlockLength(Nblock), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 45')
      DO iBlock=1,Nblock
        IF (iBlock <= Delta) THEN
          lenBlock=Hlen+1
        ELSE
          lenBlock=Hlen
        END IF
        LocalColor % BlockLength(iBlock)=lenBlock
      END DO
      LocalColor % maxBlockLength = maxBlockLength

      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'FINISHING INIT_BLOCK_FREQDIR'
      FLUSH(STAT%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_BLK_L2U_ARRAY(LocalColor)
      USE DATAPOOL, only : MNP, MSC, MDC, LocalColorInfo, rkind, stat
      USE DATAPOOL, only : wwm_nnbr_send, wwm_nnbr_recv
      USE DATAPOOL, only : wwm_ListNbCommon_send, wwm_ListNbCommon_recv
      USE DATAPOOL, only : wwm_ListDspl_send, wwm_ListDspl_recv
      USE DATAPOOL, only : wwm_ListNeigh_send, wwm_ListNeigh_recv
      USE DATAPOOL, only : XP, YP, rtype, ierr, myrank, iplg
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      integer, allocatable :: ListNeed(:), IdxRev(:)
      integer nbNeedSend_blk, nbNeedRecv_blk
      integer idx, IP
      integer nbUpp_send, nbLow_recv, iUpp, iLow, iRank
      integer maxBlockLength, idxSend, idxRecv
      integer I, IC, nbCommon, eFirst
      integer ListFirstCommon_send(wwm_nnbr_send)
      integer ListFirstCommon_recv(wwm_nnbr_recv)
      integer, allocatable :: dspl_send(:), dspl_recv(:)
      integer istat

      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'ENTERING INIT_BLK_L2U_ARRAY'
      FLUSH(STAT%FHNDL)

      maxBlockLength=LocalColor % maxBlockLength
      ListFirstCommon_send=0
      DO I=2,wwm_nnbr_send
        ListFirstCommon_send(I)=ListFirstCommon_send(I-1)+wwm_ListNbCommon_send(I-1)
      END DO
      ListFirstCommon_recv=0
      DO I=2,wwm_nnbr_recv
        ListFirstCommon_recv(I)=ListFirstCommon_recv(I-1)+wwm_ListNbCommon_recv(I-1)
      END DO
      nbUpp_send=LocalColor % nbUpp_send
      nbLow_recv=LocalColor % nbLow_recv
      allocate(LocalColor % l2u_p2dsend_type(nbUpp_send), LocalColor % l2u_p2drecv_type(nbLow_recv), LocalColor % l2u_ListNeigh_send(nbUpp_send), LocalColor % l2u_ListNeigh_recv(nbLow_recv), ListNeed(MNP), IdxRev(MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 46')
      ListNeed=0
      IdxRev=0
      nbNeedSend_blk=0
      DO iUpp=1,nbUpp_send
        I=LocalColor % ListIdxUpper_send(iUpp)
        iRank=wwm_ListNeigh_send(I)-1
        LocalColor % l2u_ListNeigh_send(iUpp)=iRank
        nbCommon=wwm_ListNbCommon_send(I)
        eFirst=ListFirstCommon_send(I)
        DO IC=1,nbCommon
          idxSend=wwm_ListDspl_send(eFirst+IC)
          IF (ListNeed(idxSend) .eq. 0) THEN
            ListNeed(idxSend)=1
            nbNeedSend_blk=nbNeedSend_blk+1
          END IF
        END DO
      END DO
      LocalColor % nbNeedSend_blk=nbNeedSend_blk
      allocate(LocalColor % IdxSend_blk(nbNeedSend_blk), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 47')
      idx=0
      DO IP=1,MNP
        IF (ListNeed(IP) .eq. 1) THEN
          idx=idx+1
          LocalColor % IdxSend_blk(idx)=IP
          IdxRev(IP)=idx
        END IF
      END DO
      DO iUpp=1,nbUpp_send
        I=LocalColor % ListIdxUpper_send(iUpp)
        nbCommon=wwm_ListNbCommon_send(I)
        eFirst=ListFirstCommon_send(I)
        ALLOCATE(dspl_send(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 48')
        DO IC=1,nbCommon
          idxSend=wwm_ListDspl_send(eFirst+IC)
          dspl_send(IC)=(IdxRev(idxSend)-1)*maxBlockLength
        END DO
        call mpi_type_create_indexed_block(nbCommon,maxBlockLength,dspl_send,rtype,LocalColor % l2u_p2dsend_type(iUpp),ierr)
        call mpi_type_commit(LocalColor % l2u_p2dsend_type(iUpp), ierr)
        DEALLOCATE(dspl_send)
      END DO
      !
      !
      ListNeed=0
      IdxRev=0
      nbNeedRecv_blk=0
      DO iLow=1,nbLow_recv
        I=LocalColor % ListIdxLower_recv(iLow)
        iRank=wwm_ListNeigh_recv(I)-1
        LocalColor % l2u_ListNeigh_recv(iLow)=iRank
        nbCommon=wwm_ListNbCommon_recv(I)
        eFirst=ListFirstCommon_recv(I)
        DO IC=1,nbCommon
          idxRecv=wwm_ListDspl_recv(eFirst+IC)
          IF (ListNeed(idxRecv) .eq. 0) THEN
            ListNeed(idxRecv)=1
            nbNeedRecv_blk=nbNeedRecv_blk+1
          END IF
        END DO
      END DO
      LocalColor % nbNeedRecv_blk=nbNeedRecv_blk
      allocate(LocalColor % IdxRecv_blk(nbNeedRecv_blk), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 49')
      idx=0
      DO IP=1,MNP
        IF (ListNeed(IP) .eq. 1) THEN
          idx=idx+1
          LocalColor % IdxRecv_blk(idx)=IP
          IdxRev(IP)=idx
        END IF
      END DO
      DO iLow=1,nbLow_recv
        I=LocalColor % ListIdxLower_recv(iLow)
        nbCommon=wwm_ListNbCommon_recv(I)
        eFirst=ListFirstCommon_recv(I)
        ALLOCATE(dspl_recv(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 50')
        DO IC=1,nbCommon
          idxRecv=wwm_ListDspl_recv(eFirst+IC)
          dspl_recv(IC)=(IdxRev(idxRecv)-1)*maxBlockLength
        END DO
        call mpi_type_create_indexed_block(nbCommon,maxBlockLength,dspl_recv,rtype,LocalColor % l2u_p2drecv_type(iLow),ierr)
        call mpi_type_commit(LocalColor % l2u_p2drecv_type(iLow), ierr)
        DEALLOCATE(dspl_recv)
      END DO
      deallocate(ListNeed, IdxRev)

      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'FINISHING INIT_BLK_L2U_ARRAY'
      FLUSH(STAT%FHNDL)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SYMM_INIT_COLORING(LocalColor, NbBlock)
      USE DATAPOOL, only : LocalColorInfo, MNP, rkind, XP, YP, stat
      USE DATAPOOL, only : DO_SYNC_UPP_2_LOW, DO_SYNC_LOW_2_UPP, DO_SYNC_FINAL
      USE datapool, only : myrank, nproc, iplg, Graph
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      integer, intent(in) :: NbBlock
      type(Graph) :: AdjGraph
      integer :: ListColor(nproc)
      integer :: ListColorWork(nproc)
      integer istat
# ifdef DEBUG
      integer TheRes
# endif

      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'ENTERING SYMM_INIT_COLORING'
      FLUSH(STAT%FHNDL)
# ifdef DEBUG
      CALL COMPUTE_TOTAL_INDEX_SHIFT(TheRes)
      WRITE(740+myrank,*) 'Total residual shift=', TheRes
# endif
      CALL COLLECT_ALL_IA_JA
# ifdef DEBUG
      WRITE(740+myrank,*) 'After COLLECT_ALL_IA_JA'
# endif
      CALL CREATE_WWM_P2D_EXCH
# ifdef DEBUG
      WRITE(740+myrank,*) 'After CREATE_WWM_P2D_EXCH'
# endif
# ifdef DEBUG
      CALL CHECK_I5B_EXCHANGE(LocalColor)
# endif
      CALL CREATE_WWM_MAT_P2D_EXCH
# ifdef DEBUG
      WRITE(740+myrank,*) 'After CREATE_WWM_MAT_P2D_EXCH'
# endif
      CALL SYMM_GRAPH_BUILD_ADJ(AdjGraph)
      CALL BUILD_MULTICOLORING(AdjGraph, ListColor)
# ifdef DEBUG
      WRITE(740+myrank,*) 'After BUILD_MULTICOLORING'
# endif
      CALL DeallocateGraph(AdjGraph)
      ListColorWork=-ListColor
      allocate(LocalColor % ListColor(nproc), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 51')
      LocalColor % ListColor=ListColorWork
# ifdef DEBUG
      WRITE(740+myrank,*) 'Before INIT_LOW_2_UPP_ARRAYS'
# endif
      CALL INIT_LOW_2_UPP_ARRAYS(LocalColor, ListColorWork)
      !
# ifdef DEBUG
      WRITE(740+myrank,*) 'Before CALL_BLOCK_FREQDIR'
# endif
      CALL INIT_BLOCK_FREQDIR(LocalColor, NbBlock)
      !
# ifdef DEBUG
      WRITE(740+myrank,*) 'Before INIT_BLK_L2U_ARRAY'
# endif
      CALL INIT_BLK_L2U_ARRAY(LocalColor)
      !
# ifdef DEBUG
      WRITE(740+myrank,*) 'Before COLLECT_ALL_COVLOWER'
# endif
      CALL COLLECT_ALL_COVLOWER(LocalColor)
      !
# ifdef DEBUG
      WRITE(740+myrank,*) 'Before INIT_COVLOWER_ARRAY'
# endif
      CALL INIT_COVLOWER_ARRAY(LocalColor)
      !
# ifdef DEBUG
      WRITE(740+myrank,*) 'Before DETERMINE_JSTATUS_L_U'
# endif
      CALL DETERMINE_JSTATUS_L_U(LocalColor)
      !
      DO_SYNC_UPP_2_LOW=.TRUE.
      DO_SYNC_LOW_2_UPP=.TRUE.
      DO_SYNC_FINAL=.TRUE.

      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'FINISHED WITH SYMM_INIT_COLORING'
      FLUSH(STAT%FHNDL)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_LOW_2_UPP_ARRAYS(LocalColor, ListColor)
      USE DATAPOOL, only : LocalColorInfo, MNP, rkind, stat
      USE DATAPOOL, only : NNZ, IA, JA, NP_RES
      USE DATAPOOL, only : wwm_ListNbCommon_send, wwm_ListNbCommon_recv
      USE DATAPOOL, only : wwm_nnbr_send, wwm_nnbr_recv
      USE DATAPOOL, only : wwm_p2drecv_type, wwm_p2dsend_type
      USE DATAPOOL, only : wwm_ListNeigh_send, wwm_ListNeigh_recv
      USE DATAPOOL, only : wwm_ListDspl_recv
      USE datapool, only : myrank, nproc, comm, ierr, nbrrank_p
      USE datapool, only : MPI_STATUS_SIZE
      implicit none

      type(LocalColorInfo), intent(inout) :: LocalColor
      integer, intent(in) :: ListColor(nproc)
      real(rkind) :: p2d_data_send(MNP)
      real(rkind) :: CovLower(MNP), CovLower_meth2(MNP)
      real(rkind) :: SumErr
      integer eColor, fColor, I, iRank
      integer nbUpp_send, nbLow_recv
      integer iLow, iUpp
      integer IC, eFirst, nbCommon, IPloc
      integer ListFirstCommon_send(wwm_nnbr_send)
      integer ListFirstCommon_recv(wwm_nnbr_recv)
      integer istat
# ifdef DEBUG
      integer IP
# endif

      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'ENTERING INIT_LOW_2_UPP_ARRAYS'
      FLUSH(STAT%FHNDL)

      eColor=ListColor(myrank+1)
      nbUpp_send=0
      DO I=1,wwm_nnbr_send
        iRank=wwm_ListNeigh_send(I)
        fColor=ListColor(iRank)
# ifdef DEBUG
        WRITE(740+myrank,*) 'I=', I, 'iRank=', iRank, 'fColor=', fColor
        IF (fColor.eq.eColor) THEN
          call wwm_abort('Major error in the code')
        END IF
# endif
        IF (fColor.gt.eColor) THEN
          nbUpp_send=nbUpp_send + 1
        ENDIF
      END DO
      nbLow_recv=0
      DO I=1,wwm_nnbr_recv
        iRank=wwm_ListNeigh_recv(I)
        fColor=ListColor(iRank)
# ifdef DEBUG
        IF (fColor.eq.eColor) THEN
          call wwm_abort('Major error in the code')
        END IF
# endif
        IF (fColor.lt.eColor) THEN
          nbLow_recv=nbLow_recv + 1
        ENDIF
      END DO
      ListFirstCommon_send=0
      DO I=2,wwm_nnbr_send
        ListFirstCommon_send(I)=ListFirstCommon_send(I-1)+wwm_ListNbCommon_send(I-1)
      END DO
      ListFirstCommon_recv=0
      DO I=2,wwm_nnbr_recv
        ListFirstCommon_recv(I)=ListFirstCommon_recv(I-1)+wwm_ListNbCommon_recv(I-1)
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'SIC: nbLow_recv=', nbLow_recv, ' nbUpp_send=', nbUpp_send
# endif
      LocalColor % nbUpp_send=nbUpp_send
      LocalColor % nbLow_recv=nbLow_recv
      allocate(LocalColor % ListIdxUpper_send(nbUpp_send), LocalColor % ListIdxLower_recv(nbLow_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 52')
      iUpp=0
      DO I=1,wwm_nnbr_send
        iRank=wwm_ListNeigh_send(I)
        fColor=ListColor(iRank)
        IF (fColor.gt.eColor) THEN
          iUpp=iUpp + 1
          LocalColor % ListIdxUpper_send(iUpp)=I
        ENDIF
      END DO
      iLow=0
      DO I=1,wwm_nnbr_recv
        iRank=wwm_ListNeigh_recv(I)
        fColor=ListColor(iRank)
        IF (fColor.lt.eColor) THEN
          iLow=iLow + 1
          LocalColor % ListIdxLower_recv(iLow)=I
        ENDIF
      END DO
      allocate(LocalColor % Upp_s_rq(nbUpp_send), LocalColor % Upp_s_stat(MPI_STATUS_SIZE, nbUpp_send), LocalColor % Low_r_rq(nbLow_recv), LocalColor % Low_r_stat(MPI_STATUS_SIZE, nbLow_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 53')
      p2d_data_send=0
      CovLower=1
      CovLower_meth2=1
      DO iUpp=1,nbUpp_send
        I=LocalColor % ListIdxUpper_send(iUpp)
        iRank=wwm_ListNeigh_send(i)
        call mpi_isend(p2d_data_send,1,wwm_p2dsend_type(i),iRank-1,13,comm,LocalColor % Upp_s_rq(iUpp),ierr)
      END DO
      DO iLow=1,nbLow_recv
        I=LocalColor % ListIdxLower_recv(iLow)
        iRank=wwm_ListNeigh_recv(i)
        nbCommon=wwm_ListNbCommon_recv(I)
        eFirst=ListFirstCommon_recv(I)
        DO IC=1,nbCommon
          IPloc=wwm_ListDspl_recv(eFirst+IC)
          CovLower_meth2(IPloc)=0
        END DO
        call mpi_irecv(CovLower,1,wwm_p2drecv_type(i),iRank-1,13,comm,LocalColor % Low_r_rq(iLow),ierr)
      END DO
      IF (nbUpp_send > 0) THEN
        call mpi_waitall(nbUpp_send, LocalColor%Upp_s_rq, LocalColor%Upp_s_stat,ierr)
      END IF
      IF (nbLow_recv > 0) THEN
        call mpi_waitall(nbLow_recv, LocalColor%Low_r_rq, LocalColor%Low_r_stat,ierr)
      END IF
      allocate(LocalColor % CovLower(MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 54')
      LocalColor % CovLower=INT(CovLower)
      SumErr=sum(abs(CovLower-CovLower_meth2))
# ifdef DEBUG
      WRITE(740+myrank,*) 'SumErr(meth1/meth2) CovLower=', SumErr
      WRITE(740+myrank,*) 'MNP=', MNP, ' sum(CovLower)=', sum(CovLower)
      IF (SumErr .gt. 0) THEN
        DO IP=1,MNP
          WRITE(740+myrank,*) IP, CovLower(IP), CovLower_meth2(IP)
        END DO
      ENDIF
# endif
      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'FINISHED WITH INIT_LOW_2_UPP_ARRAYS'
      FLUSH(STAT%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_COVLOWER_ARRAY(LocalColor)
      USE DATAPOOL
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      integer :: ListFirst(nproc)
      integer :: ListCommon_recv(nproc)
      integer :: ListCommon_send(nproc)
      integer :: ListMapped0(np_global)
      integer :: ListMapped1(np_global)
      integer :: ListMapped0_B(np_global)
      integer :: ListMapped1_B(np_global)
      integer, allocatable :: dspl_send(:), dspl_recv(:)
      integer IP, IP_glob, iProc, WeMatch, MNPloc, idx, NP_RESloc
      integer iNeigh, IPmap, nbCommon
      integer nbCommon_send, nbCommon_recv, idx_send, idx_recv
      integer u2l_nnbr_send, u2l_nnbr_recv
      integer sync_nnbr_send, sync_nnbr_recv
      integer eCov, eColor, fColor
      integer maxBlockLength
      integer nbMap0, nbMap1
      integer DoOper
      integer, allocatable :: ListNeedSend(:), ListNeedRecv(:), IdxRev(:)
      integer nbNeedSend_u2l, nbNeedRecv_u2l
      integer lenMNP

      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'ENTERING INIT_COVLOWER_ARRAY'
      FLUSH(STAT%FHNDL)

      ListFirst=0
      DO iProc=2,nproc
        ListFirst(iProc)=ListFirst(iProc-1) + ListMNP(iProc-1)
      END DO
      ListMapped0=0
      ListMapped1=0
      nbMap0=0
      nbMap1=0
      DO IP=1,MNP
        IP_glob=iplg(IP)
        eCov=LocalColor % CovLower(IP)
        IF (eCov == 0) THEN
          ListMapped0(IP_glob)=IP
          nbMap0=nbMap0+1
        END IF
        IF (eCov == 1) THEN
          ListMapped1(IP_glob)=IP
          nbMap1=nbMap1+1
        END IF
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'nbMap0=', nbMap0, ' nbMap1=', nbMap1
# endif
      !
      ! First the Upper to lower (u2l) block arrays 
      !
      u2l_nnbr_send=0
      u2l_nnbr_recv=0
      ListCommon_send=0
      ListCommon_recv=0
      eColor=LocalColor % ListColor(myrank+1)
      allocate(ListNeedRecv(MNP), ListNeedSend(MNP), IdxRev(MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 55')
      ListNeedRecv=0
      ListNeedSend=0
      IdxRev=0
      nbNeedSend_u2l=0
      nbNeedRecv_u2l=0
# ifdef DEBUG
      WRITE(740+myrank,*) 'U2L eColor=', eColor
# endif
      DO iNeigh=1,wwm_nnbr
        iProc=wwm_ListNeigh(iNeigh)
        fColor=LocalColor % ListColor(iProc)
# ifdef DEBUG
        WRITE(740+myrank,*) 'U2L iNeigh=', iNeigh, ' iProc=', iProc, 'fColor=', fColor
# endif
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        IF (fColor .ge. eColor) THEN
          nbCommon_recv=0
          DO IP=1,NP_RESloc
            IP_glob=ListIPLG(IP+ListFirst(iProc))
            eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
            IF (eCov .eq. 1) THEN
              IPmap=ListMapped0(IP_glob)
              WeMatch=0
              IF (IPmap .gt. 0) THEN
                WeMatch=1
              ELSE
                IPmap=ListMapped1(IP_glob)
                IF ((IPmap .gt. 0).and.(IPmap .gt. NP_RES)) THEN
                  WeMatch=1
                END IF
              END IF
              IF (WeMatch .eq. 1) THEN
                nbCommon_recv=nbCommon_recv+1
                IF (ListNeedRecv(IPmap) .eq. 0) THEN
                  nbNeedRecv_u2l=nbNeedRecv_u2l+1
                  ListNeedRecv(IPmap)=1
                END IF
              END IF
            END IF
          END DO
          IF (nbCommon_recv .gt. 0) THEN
            u2l_nnbr_recv=u2l_nnbr_recv+1
            ListCommon_recv(iProc)=nbCommon_recv
          END IF
# ifdef DEBUG
          WRITE(740+myrank,*) '   U2L nbCommon_recv=', nbCommon_recv
# endif
        END IF
        IF (fColor .le. eColor) THEN
          nbCommon_send=0
          DO IP=1,MNPloc
            IP_glob=ListIPLG(IP+ListFirst(iProc))
            eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
            IPmap=ListMapped1(IP_glob)
            IF ((IPmap .gt. 0).and.(IPmap .le. NP_RES)) THEN
              IF ((eCov .eq. 0).or.(IP.gt.NP_RESloc)) THEN
                nbCommon_send=nbCommon_send+1
                IF (ListNeedSend(IPmap) .eq. 0) THEN
                  nbNeedSend_u2l=nbNeedSend_u2l+1
                  ListNeedSend(IPmap)=1
                END IF
              END IF
            END IF
          END DO
          IF (nbCommon_send .gt. 0) THEN
            u2l_nnbr_send=u2l_nnbr_send+1
            ListCommon_send(iProc)=nbCommon_send
          END IF
# ifdef DEBUG
          WRITE(740+myrank,*) '   U2L nbCommon_send=', nbCommon_send
# endif
        END IF
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_P2D: u2l_nnbr_send=', u2l_nnbr_send
      WRITE(740+myrank,*) 'WWM_P2D: u2l_nnbr_recv=', u2l_nnbr_recv
# endif
      LocalColor % u2l_nnbr_send=u2l_nnbr_send
      LocalColor % u2l_nnbr_recv=u2l_nnbr_recv
      allocate(LocalColor % u2l_ListNbCommon_send(u2l_nnbr_send), LocalColor % u2l_ListNbCommon_recv(u2l_nnbr_recv), LocalColor % u2l_ListNeigh_send(u2l_nnbr_send), LocalColor % u2l_ListNeigh_recv(u2l_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 56')
      idx_send=0
      idx_recv=0
      DO iProc=1,nproc
        IF (ListCommon_send(iProc) .gt. 0) THEN
          idx_send=idx_send+1
          LocalColor % u2l_ListNeigh_send(idx_send)=iProc-1
          nbCommon=ListCommon_send(iProc)
          LocalColor % u2l_ListNbCommon_send(idx_send)=nbCommon
        END IF
        IF (ListCommon_recv(iProc) .gt. 0) THEN
          idx_recv=idx_recv+1
          LocalColor % u2l_ListNeigh_recv(idx_recv)=iProc-1
          nbCommon=ListCommon_recv(iProc)
          LocalColor % u2l_ListNbCommon_recv(idx_recv)=nbCommon
        END IF
      END DO
      LocalColor % nbNeedSend_u2l=nbNeedSend_u2l
      LocalColor % nbNeedRecv_u2l=nbNeedRecv_u2l
      allocate(LocalColor % IdxSend_u2l(nbNeedSend_u2l), LocalColor % IdxRecv_u2l(nbNeedRecv_u2l), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 57')
      !
      ! Now creating the u2l exchange
      ! 
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_P2D: wwm_nnbr=', wwm_nnbr
      WRITE(740+myrank,*) 'WWM_P2D: wwm_ListNeigh built'
# endif
      allocate(LocalColor % u2l_p2dsend_rqst(u2l_nnbr_send), LocalColor % u2l_p2drecv_rqst(u2l_nnbr_recv), &
     &LocalColor % u2l_p2dsend_stat(MPI_STATUS_SIZE,u2l_nnbr_send), &
     &LocalColor % u2l_p2drecv_stat(MPI_STATUS_SIZE,u2l_nnbr_recv), &
     &LocalColor % u2l_p2dsend_type(u2l_nnbr_send), LocalColor % u2l_p2drecv_type(u2l_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 58')
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_P2D: alloc done'
# endif
      maxBlockLength=LocalColor % maxBlockLength
      idx=0
      DO IP=1,MNP
        IF (ListNeedSend(IP) .eq. 1) THEN
          idx=idx+1
          LocalColor % IdxSend_u2l(idx)=IP
          IdxRev(IP)=idx
        END IF
      END DO
      DO iNeigh=1,u2l_nnbr_send
        iProc=LocalColor % u2l_ListNeigh_send(iNeigh)+1
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        ListMapped0_B=0
        ListMapped1_B=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 0) THEN
            ListMapped0_B(IP_glob)=IP
          END IF
          IF (eCov == 1) THEN
            ListMapped1_B(IP_glob)=IP
          END IF
        END DO
        nbCommon=LocalColor % u2l_ListNbCommon_send(iNeigh)
        allocate(dspl_send(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 59')
        idx=0
        DO IP=1,NP_RES
          IF (LocalColor % CovLower(IP) .eq. 1) THEN
            IP_glob=iplg(IP)
            IPmap=ListMapped0_B(IP_glob)
            WeMatch=0
            IF (IPmap .gt. 0) THEN
              WeMatch=1
            ELSE
              IPmap=ListMapped1_B(IP_glob)
              IF ((IPmap .gt. 0).and.(IPmap .gt. NP_RESloc)) THEN
                WeMatch=1
              END IF
            END IF
            IF (WeMatch .eq. 1) THEN
              idx=idx+1
              dspl_send(idx)=maxBlockLength*(IdxRev(IP)-1)
            END IF
          END IF
        END DO
        call mpi_type_create_indexed_block(nbCommon,maxBlockLength,dspl_send,rtype,LocalColor % u2l_p2dsend_type(iNeigh), ierr)
        call mpi_type_commit(LocalColor % u2l_p2dsend_type(iNeigh), ierr)
        deallocate(dspl_send)
      END DO
      IdxRev=0
      idx=0
      DO IP=1,MNP
        IF (ListNeedRecv(IP) .eq. 1) THEN
          idx=idx+1
          LocalColor % IdxRecv_u2l(idx)=IP
          IdxRev(IP)=idx
        END IF
      END DO
      DO iNeigh=1,u2l_nnbr_recv
        iProc=LocalColor % u2l_ListNeigh_recv(iNeigh)+1
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        ListMapped0_B=0
        ListMapped1_B=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 0) THEN
            ListMapped0_B(IP_glob)=IP
          END IF
          IF (eCov == 1) THEN
            ListMapped1_B(IP_glob)=IP
          END IF
        END DO
        nbCommon=LocalColor % u2l_ListNbCommon_recv(iNeigh)
        allocate(dspl_recv(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 60')
        idx=0
        DO IP=1,NP_RESloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 1) THEN
            IPmap=ListMapped0(IP_glob)
            WeMatch=0
            IF (IPmap .gt. 0) THEN
              WeMatch=1
            ELSE
              IPmap=ListMapped1(IP_glob)
              IF ((IPmap .gt. 0).and.(IPmap .gt. NP_RES)) THEN
                WeMatch=1
              END IF
            END IF
            IF (WeMatch .eq. 1) THEN
              idx=idx+1
              dspl_recv(idx)=maxBlockLength*(IdxRev(IPmap)-1)
            END IF
          END IF
        END DO
        call mpi_type_create_indexed_block(nbCommon,maxBlockLength,dspl_recv,rtype,LocalColor % u2l_p2drecv_type(iNeigh),ierr)
        call mpi_type_commit(LocalColor % u2l_p2drecv_type(iNeigh), ierr)
        deallocate(dspl_recv)
      END DO
      deallocate(ListNeedRecv, ListNeedSend, IdxRev)
      !
      ! Now the synchronization arrays
      !
      sync_nnbr_send=0
      sync_nnbr_recv=0
      ListCommon_send=0
      ListCommon_recv=0
# ifdef DEBUG
      WRITE(740+myrank,*) 'wwm_nnbr=', wwm_nnbr
      WRITE(740+myrank,*) 'sum(ListCovLower)=', sum(LocalColor % ListCovLower)
# endif
      DO iNeigh=1,wwm_nnbr
        iProc=wwm_ListNeigh(iNeigh)
        fColor=LocalColor % ListColor(iProc)
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        ListMapped0_B=0
        ListMapped1_B=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 0) THEN
            ListMapped0_B(IP_glob)=IP
          END IF
          IF (eCov == 1) THEN
            ListMapped1_B(IP_glob)=IP
          END IF
        END DO
        nbCommon_recv=0
        DO IP=1,NP_RESloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov .eq. 1) THEN
            IF (ListMapped0(IP_glob) .gt. 0) THEN
              nbCommon_recv=nbCommon_recv+1
            ELSE
              IPmap=ListMapped1(IP_glob)
              IF ((IPmap .gt. 0).and.(IPmap .gt. NP_RES)) THEN
                nbCommon_recv=nbCommon_recv+1
              ENDIF
            END IF
          END IF
        END DO
        IF (nbCommon_recv .gt. 0) THEN
          sync_nnbr_recv=sync_nnbr_recv+1
          ListCommon_recv(iProc)=nbCommon_recv
        END IF
        nbCommon_send=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IPmap=ListMapped1(IP_glob)
          IF ((IPmap .gt. 0).and.(IPmap .le. NP_RES)) THEN
            IF ((eCov .eq. 0).or.(IP.gt.NP_RESloc)) THEN
              nbCommon_send=nbCommon_send+1
            END IF
          END IF
        END DO
        IF (nbCommon_send .gt. 0) THEN
          sync_nnbr_send=sync_nnbr_send+1
          ListCommon_send(iProc)=nbCommon_send
        END IF
# ifdef DEBUG
        WRITE(740+myrank,*) '   nbCommon(send/recv)=', nbCommon_send, nbCommon_recv
# endif
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_P2D: sync_nnbr_send=', sync_nnbr_send
      WRITE(740+myrank,*) 'WWM_P2D: sync_nnbr_recv=', sync_nnbr_recv
# endif
      LocalColor % sync_nnbr_send=sync_nnbr_send
      LocalColor % sync_nnbr_recv=sync_nnbr_recv
      allocate(LocalColor % sync_ListNbCommon_send(sync_nnbr_send), LocalColor % sync_ListNbCommon_recv(sync_nnbr_recv), LocalColor % sync_ListNeigh_send(sync_nnbr_send), LocalColor % sync_ListNeigh_recv(sync_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 61')
      idx_send=0
      idx_recv=0
      DO iProc=1,nproc
        IF (ListCommon_send(iProc) .gt. 0) THEN
          idx_send=idx_send+1
          LocalColor % sync_ListNeigh_send(idx_send)=iProc-1
          nbCommon=ListCommon_send(iProc)
          LocalColor % sync_ListNbCommon_send(idx_send)=nbCommon
        END IF
        IF (ListCommon_recv(iProc) .gt. 0) THEN
          idx_recv=idx_recv+1
          LocalColor % sync_ListNeigh_recv(idx_recv)=iProc-1
          nbCommon=ListCommon_recv(iProc)
          LocalColor % sync_ListNbCommon_recv(idx_recv)=nbCommon
        END IF
      END DO
      !
      ! Now creating the sync exchange
      !
      allocate(LocalColor % sync_p2dsend_rqst(sync_nnbr_send), &
     &LocalColor % sync_p2drecv_rqst(sync_nnbr_recv), &
     &LocalColor % sync_p2dsend_stat(MPI_STATUS_SIZE,sync_nnbr_send), &
     &LocalColor % sync_p2drecv_stat(MPI_STATUS_SIZE,sync_nnbr_recv), &
     &LocalColor % sync_p2dsend_type(sync_nnbr_send), &
     &LocalColor % sync_p2drecv_type(sync_nnbr_recv), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 62')
# ifdef DEBUG
      WRITE(740+myrank,*) 'SYNC sync_nnbr_send=', sync_nnbr_send
# endif
      DO iNeigh=1,sync_nnbr_send
        iProc=LocalColor % sync_ListNeigh_send(iNeigh)+1
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        ListMapped0_B=0
        ListMapped1_B=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 0) THEN
            ListMapped0_B(IP_glob)=IP
          END IF
          IF (eCov == 1) THEN
            ListMapped1_B(IP_glob)=IP
          END IF
        END DO
        nbCommon=LocalColor % sync_ListNbCommon_send(iNeigh)
# ifdef DEBUG
        WRITE(740+myrank,*) '   SYNC iNeigh=', iNeigh, ' nbCommon=', nbCommon
# endif
        allocate(dspl_send(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 63')
        idx=0
        DO IP=1,NP_RES
          IF (LocalColor % CovLower(IP) .eq. 1) THEN
            DoOper=0
            IP_glob=iplg(IP)
            IPmap=ListMapped0_B(IP_glob)
            IF (IPmap .gt. 0) THEN
              DoOper=1
            ELSE
              IPmap=ListMapped1_B(IP_glob)
              IF ((IPmap .gt. 0).and.(IPmap.gt.NP_RESloc)) THEN
                DoOper=1
              END IF
            END IF
            IF (DoOper == 1) THEN
              idx=idx+1
              dspl_send(idx)=MSC*MDC*(IP-1)
            END IF
          END IF
        END DO
        call mpi_type_create_indexed_block(nbCommon,MSC*MDC,dspl_send,rtype,LocalColor % sync_p2dsend_type(iNeigh), ierr)
        call mpi_type_commit(LocalColor % sync_p2dsend_type(iNeigh), ierr)
        deallocate(dspl_send)
      END DO
      DO iNeigh=1,sync_nnbr_recv
        iProc=LocalColor % sync_ListNeigh_recv(iNeigh)+1
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        ListMapped0_B=0
        ListMapped1_B=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 0) THEN
            ListMapped0_B(IP_glob)=IP
          END IF
          IF (eCov == 1) THEN
            ListMapped1_B(IP_glob)=IP
          END IF
        END DO
        nbCommon=LocalColor%sync_ListNbCommon_recv(iNeigh)
        allocate(dspl_recv(nbCommon), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 64')
        idx=0
        DO IP=1,NP_RESloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 1) THEN
            IPmap=ListMapped0(IP_glob)
            DoOper=0
            IF (IPmap .gt. 0) THEN
              DoOper=1
            ELSE
              IPmap=ListMapped1(IP_glob)
              IF ((IPmap .gt. 0).and.(IPmap.gt.NP_RES)) THEN
                DoOper=1
              END IF
            END IF
            IF (DoOper == 1) THEN
              idx=idx+1
              dspl_recv(idx)=MSC*MDC*(IPmap-1)
            END IF
          END IF
        END DO
        call mpi_type_create_indexed_block(nbCommon,MSC*MDC,dspl_recv,rtype,LocalColor % sync_p2drecv_type(iNeigh),ierr)
        call mpi_type_commit(LocalColor % sync_p2drecv_type(iNeigh), ierr)
        deallocate(dspl_recv)
      END DO
      lenMNP=0
      IF (LocalColor % nbNeedSend_blk > lenMNP) THEN
        lenMNP=LocalColor % nbNeedSend_blk
      END IF
      IF (LocalColor % nbNeedRecv_blk > lenMNP) THEN
        lenMNP=LocalColor % nbNeedRecv_blk
      END IF
      IF (LocalColor % nbNeedSend_u2l > lenMNP) THEN
        lenMNP=LocalColor % nbNeedSend_u2l
      END IF
      IF (LocalColor % nbNeedRecv_u2l > lenMNP) THEN
        lenMNP=LocalColor % nbNeedRecv_u2l
      END IF
      allocate(LocalColor % ACexch(maxBlockLength, lenMNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 65')

      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'FINISHED WITH INIT_COVLOWER_ARRAY'
      FLUSH(STAT%FHNDL)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DETERMINE_JSTATUS_L_U(LocalColor)
      USE DATAPOOL, only : LocalColorInfo, MNP, rkind
      USE DATAPOOL, only : NNZ, IA, JA, NP_RES, I_DIAG
      USE DATAPOOL, only : wwm_nnbr_send, wwm_nnbr_recv
      USE datapool, only : myrank, nproc, comm, ierr, nbrrank_p
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      integer Jstatus_L(NNZ), Jstatus_U(NNZ)
      integer IP, J, JP, DoOper
      integer istat
# ifdef REORDER_ASPAR_PC
      integer nb, idx
# endif
      Jstatus_L=0
      Jstatus_U=0
      DO IP=1,NP_RES
        IF (LocalColor%CovLower(IP) == 1) THEN
          DO J=IA(IP),IA(IP+1)-1
            JP=JA(J)
            IF (JP /= IP) THEN
              IF (JP .lt. IP) THEN
                DoOper=1
              ELSE
                IF (LocalColor % CovLower(JP) == 0) THEN
                  DoOper=1
                ELSE
                  DoOper=0
                END IF
              END IF
            ELSE
              DoOper=0
            END IF
            Jstatus_L(J)=DoOper
          END DO
          DO J=IA(IP),IA(IP+1)-1
            JP=JA(J)
            IF (JP /= IP) THEN
              IF (JP .lt. IP) THEN
                DoOper=0
              ELSE
                IF (LocalColor % CovLower(JP) == 0) THEN
                  DoOper=0
                ELSE
                  DoOper=1
                END IF
              END IF
            ELSE
              DoOper=0
            END IF
            Jstatus_U(J)=DoOper
          END DO
        END IF
      END DO
# ifdef DEBUG
      DO IP=1,NP_RES
        DO J=IA(IP),IA(IP+1)-1
          JP=JA(J)
          IF ((Jstatus_L(J).eq.1).and.(Jstatus_U(J).eq.1)) THEN
            WRITE(myrank+919,*) 'MNP=', MNP, 'NP_RES=', NP_RES
            WRITE(myrank+919,*) 'IP=', IP, ' JP=', JP
            WRITE(myrank+919,*) 'IPcovLower=', LocalColor%CovLower(IP)
            WRITE(myrank+919,*) 'JPcovLower=', LocalColor%CovLower(JP)
            WRITE(myrank+919,*) 'We have major error'
            CALL WWM_ABORT('Please panic and debug')
          END IF
        END DO
      END DO
      WRITE(740+myrank,*) 'sum(Jstatus_L)=', sum(Jstatus_L)
      WRITE(740+myrank,*) 'sum(Jstatus_U)=', sum(Jstatus_U)
# endif
# ifdef REORDER_ASPAR_PC
      allocate(LocalColor % IA_L(NP_RES+1), LocalColor % IA_U(NP_RES+1), LocalColor % JA_LU(NNZ), LocalColor % JmapR(NNZ), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('determine_jstatus_L_U, error 1')
      LocalColor % JmapR=-1
      LocalColor % IA_L(1)=1
      idx=0
      DO IP=1,NP_RES
        nb=0
        DO J=IA(IP),IA(IP+1)-1
          IF (Jstatus_L(J).eq.1) THEN
            JP=JA(J)
            nb=nb+1
            idx=idx+1
            LocalColor % JmapR(J)=idx
            LocalColor % JA_LU(idx)=JP
          END IF
        END DO
        LocalColor % IA_L(IP+1)=LocalColor % IA_L(IP)+nb
      END DO
      LocalColor % IA_U(1)=LocalColor % IA_L(NP_RES+1)
      DO IP=1,NP_RES
        nb=0
        DO J=IA(IP),IA(IP+1)-1
          IF (Jstatus_U(J).eq.1) THEN
            JP=JA(J)
            nb=nb+1
            idx=idx+1
            LocalColor % JmapR(J)=idx
            LocalColor % JA_LU(idx)=JP
          END IF
        END DO
        nb=nb+1
        idx=idx+1
        J=I_DIAG(IP)
        LocalColor % JmapR(J)=idx
        LocalColor % JA_LU(idx)=JP
        LocalColor % IA_U(IP+1)=LocalColor % IA_U(IP)+nb
      END DO
# endif
      allocate(LocalColor % Jstatus_L(NNZ), LocalColor % Jstatus_U(NNZ), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 66')
      LocalColor % Jstatus_L=Jstatus_L
      LocalColor % Jstatus_U=Jstatus_U
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_RECV_ASPAR_PC(LocalColor, SolDat)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : MNP, MSC, MDC, rkind
      USE datapool, only : ierr, comm, rtype, istatus, nbrrank_p
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      real(rkind), allocatable :: ASPAR_rs(:)
      integer idx, iNNZ, jNNZ, IS, ID, NNZ_l, siz
      integer iProc, i, iRank
      integer istat
      DO iProc=1,LocalColor % nbLow_send
        i=LocalColor % ListIdxUpper_send(iProc)
        iRank=nbrrank_p(i)
        NNZ_l=LocalColor % NNZ_len_r(iProc)
        siz=NNZ_l*MSC*MDC
        allocate(ASPAR_rs(NNZ_l*MSC*MDC), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 67')
        CALL mpi_recv(ASPAR_rs,siz,rtype,iRank,45,comm,istatus,ierr)
        idx=0
        DO iNNZ=1,NNZ_l
          jNNZ=LocalColor % NNZ_index_r(iProc,iNNZ)
          DO IS=1,MSC
            DO ID=1,MDC
              idx=idx+1
              SolDat % ASPAR_pc(jNNZ,IS,ID)=ASPAR_rs(idx)
            END DO
          END DO
        END DO
        deallocate(ASPAR_rs)
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_SEND_ASPAR_PC(LocalColor, SolDat)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : MSC, MDC, rkind
      USE datapool, only : ierr, comm, rtype, nbrrank_p
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      real(rkind), allocatable :: ASPAR_rs(:)
      integer idx, iNNZ, jNNZ, IS, ID, NNZ_l, siz
      integer iProc, iRank, i
      integer istat
      DO iProc=1,LocalColor % nbUpp_send
        i=LocalColor % ListIdxUpper_send(iProc)
        iRank=nbrrank_p(i)
        NNZ_l=LocalColor % NNZ_len_s(iProc)
        siz=NNZ_l*MSC*MDC
        allocate(ASPAR_rs(NNZ_l*MSC*MDC), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 68')
        idx=0
        DO iNNZ=1,NNZ_l
          jNNZ=LocalColor % NNZ_index_s(iProc,iNNZ)
          DO IS=1,MSC
            DO ID=1,MDC
              idx=idx+1
              ASPAR_rs(idx)=SolDat % ASPAR_pc(jNNZ,IS,ID)
            END DO
          END DO
        END DO
        CALL mpi_isend(ASPAR_rs,siz,rtype,iRank,45,comm,LocalColor%Upp_s_rq(iProc),ierr)
        deallocate(ASPAR_rs)
      END DO
      IF (LocalColor % nbUpp_send > 0) THEN
        call mpi_waitall(LocalColor %nbUpp_send, LocalColor % Upp_s_rq, LocalColor % Upp_s_stat,ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_CREATE_PRECOND_ILU0(LocalColor, SolDat)
      USE DATAPOOL, only : MNP, NP_RES, LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : MSC, MDC, IA, JA, I_DIAG, rkind
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      integer IP, JP, J, JP2, J2, J_FOUND, IS, ID
      integer :: ListJ(MNP)
      real(rkind) tl
      SolDat%ASPAR_pc=SolDat%ASPAR_block
      CALL I5_RECV_ASPAR_PC(LocalColor, SolDat)
      DO IP=1,NP_RES
        IF (LocalColor % CovLower(IP) == 1) THEN
          DO J=IA(IP),IA(IP+1)-1
            JP=JA(J)
            ListJ(JP)=J
          END DO
          DO J=IA(IP),I_DIAG(IP)-1
            JP=JA(J)
            DO IS=1,MSC
              DO ID=1,MDC
                tl=SolDat%ASPAR_pc(J,IS,ID)*SolDat%ASPAR_pc(I_DIAG(JP),IS,ID)
                DO J2=IA(JP),IA(JP+1)-1
                  JP2=JA(J2)
                  J_FOUND=ListJ(JP2)
                  IF (J_FOUND.gt.0) THEN ! Here is ILU0 approximation
                    SolDat%ASPAR_pc(J_FOUND,IS,ID)=SolDat%ASPAR_pc(J_FOUND,IS,ID) - tl*SolDat%ASPAR_pc(J2,IS,ID)
                  END IF
                END DO
                SolDat%ASPAR_pc(J,IS,ID)=tl
              END DO
            END DO
          END DO
          DO J=IA(IP),IA(IP+1)-1
            JP=JA(J)
            ListJ(JP)=0
          END DO
          J=I_DIAG(IP)
          DO IS=1,MSC
            DO ID=1,MDC
              SolDat%ASPAR_pc(J,IS,ID)=1.0_rkind/SolDat%ASPAR_pc(J,IS,ID)
            END DO
          END DO
        END IF
      END DO
      CALL I5_SEND_ASPAR_PC(LocalColor, SolDat)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_CREATE_PRECOND_ILU0(LocalColor, SolDat)
      USE DATAPOOL, only : MNP, NP_RES, LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : MSC, MDC, IA, JA, I_DIAG, rkind
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      CALL WWM_ABORT('Please program it')
      END SUBROUTINE
!**********************************************************************
!* We assign the values only for CovLower(IP)=1                       *
!* We could with some effort assign values for all with some effort   *
!* but the values would not be used                                   *
!**********************************************************************
# ifndef SOR_DIRECT
      SUBROUTINE I5B_CREATE_PRECOND_SOR(LocalColor, SolDat)
      USE DATAPOOL, only : MNP, MDC, IA, JA, I_DIAG, NP_RES, rkind, ONE
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      USE datapool, only : exchange_p4d_wwm
#  ifdef DEBUG
      USE datapool, only : myrank
#  endif
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      integer IP, ID, IS, JP, JR, J1, J, IPglob, JPglob
      real(rkind) eVal, Lerror
#  if defined DEBUG
      real(rkind) :: eSum, eSumB
#  endif
      DO IP=1,NP_RES
        J=I_DIAG(IP)
        SolDat%AC4(:,:,IP)=ONE/SolDat % ASPAR_block(:,:,J)
      END DO
!#  ifdef NO_SELFE_EXCH
!      CALL I5B_EXCHANGE_P4D_WWM(LocalColor, SolDat%AC4)
!#  else
!      CALL I5B_TOTAL_COHERENCY_ERROR_NPRES(SolDat%AC4, Lerror)
!      Print *, 'AC4_1 NP_RES cohenrency error=', Lerror
      CALL EXCHANGE_P4D_WWM(SolDat%AC4)
!#  endif
      DO IP=1,NP_RES
        IF (LocalColor%CovLower(IP) == 1) THEN
#  if defined REORDER_ASPAR_PC
          DO J=IA(IP),IA(IP+1)-1
            IF (LocalColor%Jstatus_L(J) == 1) THEN
              JP=JA(J)
              JR=LocalColor%JmapR(J)
              SolDat % ASPAR_pc(:,:,JR)=SolDat % ASPAR_block(:,:,J)*SolDat%AC4(:,:,JP)
            END IF
          ENDDO
          J=LocalColor% IA_U(IP+1)-1
          SolDat % ASPAR_pc(:,:,J)=SolDat%AC4(:,:,IP)
          DO J=IA(IP),IA(IP+1)-1
            IF (LocalColor%Jstatus_U(J) == 1) THEN
              JP=JA(J)
              JR=LocalColor%JmapR(J)
              SolDat % ASPAR_pc(:,:,JR)=SolDat % ASPAR_block(:,:,J)
            END IF
          END DO
#  else
          DO J=IA(IP),IA(IP+1)-1
            IF (LocalColor%Jstatus_L(J) == 1) THEN
              JP=JA(J)
              SolDat % ASPAR_pc(:,:,J)=SolDat % ASPAR_block(:,:,J)*SolDat%AC4(:,:,JP)
            END IF
          ENDDO
          J=I_DIAG(IP)
          SolDat % ASPAR_pc(:,:,J)=SolDat%AC4(:,:,IP)
          DO J=IA(IP),IA(IP+1)-1
            IF (LocalColor%Jstatus_U(J) == 1) THEN
              JP=JA(J)
              SolDat % ASPAR_pc(:,:,J)=SolDat % ASPAR_block(:,:,J)
            END IF
          END DO
#  endif
        END IF
      END DO
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_CREATE_PRECOND(LocalColor, SolDat, TheMethod)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData, I_DIAG, NP_RES, rkind
      USE datapool, only : exchange_p4d_wwm, myrank
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      integer, intent(in) :: TheMethod
      integer IP, J
!      real(rkind) :: Lerror
# ifdef SOR_DIRECT
      IF (TheMethod.ne.1) THEN
        CALL WWM_ABORT('With SOR_DIRECT, only SOR is possible')
      END IF
      DO IP=1,NP_RES
        J=I_DIAG(IP)
        SolDat%AC4(:,:,IP)=SolDat % ASPAR_block(:,:,J)
      END DO
!#  ifdef NO_SELFE_EXCH
!      CALL I5B_EXCHANGE_P4D_WWM(LocalColor, SolDat%AC4)
!#  else
!      CALL I5B_TOTAL_COHERENCY_ERROR_NPRES(SolDat%AC4, Lerror)
!      Print *, 'AC4_2 NP_RES cohenrency error=', Lerror
      CALL EXCHANGE_P4D_WWM(SolDat%AC4)
!#  endif
      DO IP=1,NP_RES
        J=I_DIAG(IP)
        SolDat % ASPAR_block(:,:,J)=SolDat%AC4(:,:,IP)
      END DO
!check this write ...
      !WRITE(myrank+640,*) 'MIN ASPAR_block=', minval(SolDat % ASPAR_block)
# else
      IF (TheMethod == 1) THEN ! SOR 
        CALL I5B_CREATE_PRECOND_SOR(LocalColor, SolDat)
      ELSE IF (TheMethod == 2) THEN ! ILU0
        CALL I5B_CREATE_PRECOND_ILU0(LocalColor, SolDat)
      ELSE
        CALL WWM_ABORT('Wrong choice of preconditioner')
      ENDIF
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_EXCHANGE_P3_LOW_2_UPP_Send(LocalColor, AC, iBlock)
      USE DATAPOOL, only : LocalColorInfo, MNP, MSC, MDC, rkind
      USE datapool, only : comm, ierr, myrank
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      REAL(rkind), intent(in) :: AC(MSC, MDC, MNP)
      INTEGER, intent(in) :: iBlock
      integer iUpp, iRank, idx, lenBlock, maxBlockLength, IS, ID, IP, nbUpp_send
      integer idxIP
      lenBlock=LocalColor % BlockLength(iBlock)
      maxBlockLength=LocalColor % maxBlockLength
      DO idxIP=1,LocalColor % nbNeedSend_blk
        IP = LocalColor % IdxSend_blk(idxIP)
        DO idx=1,lenBlock
          IS=LocalColor % ISindex(iBlock, idx)
          ID=LocalColor % IDindex(iBlock, idx)
          LocalColor % ACexch(idx,idxIP)=AC(IS,ID,IP)
        END DO
      END DO
      nbUpp_send=LocalColor % nbUpp_send
      DO iUpp=1,nbUpp_send
        iRank=LocalColor % l2u_ListNeigh_send(iUpp)
        CALL mpi_isend(LocalColor % ACexch, 1, LocalColor % l2u_p2dsend_type(iUpp), iRank, 7, comm, LocalColor%Upp_s_rq(iUpp), ierr)
      END DO
      IF (nbUpp_send > 0) THEN
        call mpi_waitall(nbUpp_send, LocalColor % Upp_s_rq, LocalColor % Upp_s_stat,ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_EXCHANGE_P3_LOW_2_UPP_Recv(LocalColor, AC, iBlock)
      USE DATAPOOL, only : LocalColorInfo, MNP, MSC, MDC, rkind
      USE datapool, only : comm, ierr, myrank
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      REAL(rkind), intent(inout) :: AC(MSC, MDC, MNP)
      INTEGER, intent(in) :: iBlock
      integer iProc, iRank, idx, lenBlock, IS, ID, IP, nbLow_recv
      integer idxIP
      lenBlock=LocalColor % BlockLength(iBlock)
      nbLow_recv=LocalColor % nbLow_recv
      DO iProc=1,nbLow_recv
        iRank=LocalColor % l2u_ListNeigh_recv(iProc)
        call mpi_irecv(LocalColor % ACexch,1,LocalColor % l2u_p2drecv_type(iProc),iRank,7,comm,LocalColor % Low_r_rq(iProc),ierr)
      END DO
      IF (nbLow_recv > 0) THEN
        call mpi_waitall(nbLow_recv, LocalColor%Low_r_rq, LocalColor%Low_r_stat,ierr)
      END IF
      DO idxIP=1,LocalColor % nbNeedRecv_blk
        IP = LocalColor % IdxRecv_blk(idxIP)
        DO idx=1,lenBlock
          IS=LocalColor % ISindex(iBlock, idx)
          ID=LocalColor % IDindex(iBlock, idx)
          AC(IS,ID,IP) = LocalColor % ACexch(idx,idxIP)
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_EXCHANGE_P3_UPP_2_LOW_Send(LocalColor, AC, iBlock)
      USE DATAPOOL, only : LocalColorInfo, MNP, MSC, MDC, rkind
      USE datapool, only : comm, ierr, myrank
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      REAL(rkind), intent(in) :: AC(MSC, MDC, MNP)
      INTEGER, intent(in) :: iBlock
      integer iProc, iRank, idx, lenBlock, IS, ID, nbLow_send, IP
      integer idxIP
      lenBlock=LocalColor % BlockLength(iBlock)
      DO idxIP=1,LocalColor % nbNeedSend_u2l
        IP = LocalColor % IdxSend_u2l(idxIP)
        DO idx=1,lenBlock
          IS=LocalColor % ISindex(iBlock, idx)
          ID=LocalColor % IDindex(iBlock, idx)
          LocalColor % ACexch(idx,idxIP)=AC(IS,ID,IP)
        END DO
      END DO
      nbLow_send=LocalColor % u2l_nnbr_send
      DO iProc=1,nbLow_send
        iRank=LocalColor % u2l_ListNeigh_send(iProc)
        call mpi_isend(LocalColor % ACexch,1,LocalColor%u2l_p2dsend_type(iProc),iRank,1151,comm,LocalColor%u2l_p2dsend_rqst(iProc),ierr)
      END DO
      IF (nbLow_send > 0) THEN
        call mpi_waitall(nbLow_send, LocalColor%u2l_p2dsend_rqst, LocalColor%u2l_p2dsend_stat,ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_EXCHANGE_P3_UPP_2_LOW_Recv(LocalColor, AC, iBlock)
      USE DATAPOOL, only : LocalColorInfo, MNP, MSC, MDC, rkind
      USE datapool, only : comm, ierr, myrank
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      REAL(rkind), intent(inout) :: AC(MSC, MDC, MNP)
      INTEGER, intent(in) :: iBlock
      integer iProc, iRank, idx, lenBlock
      integer nbUpp_recv, IS, ID, IP
      integer idxIP
      lenBlock=LocalColor % BlockLength(iBlock)
      nbUpp_recv=LocalColor % u2l_nnbr_recv
      DO iProc=1,nbUpp_recv
        iRank=LocalColor % u2l_ListNeigh_recv(iProc)
        call mpi_irecv(LocalColor % ACexch,1,LocalColor%u2l_p2drecv_type(iProc),iRank,1151,comm,LocalColor % u2l_p2drecv_rqst(iProc),ierr)
      END DO
      IF (nbUpp_recv > 0) THEN
        call mpi_waitall(nbUpp_recv, LocalColor%u2l_p2drecv_rqst, LocalColor%u2l_p2drecv_stat,ierr)
      END IF
      DO idxIP=1,LocalColor % nbNeedRecv_u2l
        IP = LocalColor % IdxRecv_u2l(idxIP)
        DO idx=1,lenBlock
          IS=LocalColor % ISindex(iBlock, idx)
          ID=LocalColor % IDindex(iBlock, idx)
          AC(IS,ID,IP) = LocalColor % ACexch(idx,idxIP)
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
# ifdef DEBUG
      SUBROUTINE COMPUTE_TOTAL_INDEX_SHIFT(TheRes)
      USE DATAPOOL, only : NP_RES, IA, JA
      implicit none
      integer, intent(out) :: TheRes
      integer :: IP, JP, J
      TheRes=0
      DO IP=1,NP_RES
        DO J=IA(IP),IA(IP+1)-1
          JP=JA(J)
          TheRes=TheRes+abs(IP-JP)
        END DO
      END DO
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_PARTIAL_SOLVE_L(LocalColor, SolDat, iBlock, ACret)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData, I_DIAG, ONE
      USE DATAPOOL, only : IA, JA, MSC, MDC, MNP, rkind, NP_RES, THR, THR8
      USE datapool, only : myrank
# ifdef DEBUG
      USE datapool, only : iplg
# endif
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(in) :: SolDat
      integer, intent(in) :: iBlock
      real(rkind), intent(inout) :: ACret(MSC,MDC,MNP)
      real(rkind) :: eCoeff
!      real(rkind) :: hVal
      integer IP, idx, ID, IS, J
      integer lenBlock, JP, Jb
!      real(rkind) :: ACtest(MSC, MDC,MNP)
      lenBlock=LocalColor % BlockLength(iBlock)
!      DO IP=1,NP_RES
!        J=I_DIAG(IP)
!        ACtest(:,:,IP)=ONE/SolDat % ASPAR_block(:,:,J)
!      END DO
      DO IP=1,NP_RES
        IF (LocalColor % CovLower(IP) == 1) THEN
# if defined REORDER_ASPAR_PC
          DO J=LocalColor % IA_L(IP),LocalColor % IA_L(IP+1)-1
            JP=LocalColor % JA_LU(J)
            DO idx=1,lenBlock
              IS=LocalColor % ISindex(iBlock, idx)
              ID=LocalColor % IDindex(iBlock, idx)
              eCoeff=SolDat % ASPAR_pc(IS,ID,J)
              ACret(IS,ID,IP)=ACret(IS,ID,IP) - eCoeff*ACret(IS,ID,JP)
            END DO
          END DO
# elif defined SOR_DIRECT
          DO J=IA(IP),IA(IP+1)-1
            IF (LocalColor % Jstatus_L(J) .eq. 1) THEN
              JP=JA(J)
              Jb=I_DIAG(JP)
              DO idx=1,lenBlock
                IS=LocalColor % ISindex(iBlock, idx)
                ID=LocalColor % IDindex(iBlock, idx)
                eCoeff=SolDat % ASPAR_block(IS,ID,J)/SolDat % ASPAR_block(IS,ID,Jb)
                ACret(IS,ID,IP)=ACret(IS,ID,IP) - eCoeff*ACret(IS,ID,JP)
              END DO
            END IF
          END DO
# else
          DO J=IA(IP),IA(IP+1)-1
            IF (LocalColor % Jstatus_L(J) .eq. 1) THEN
              JP=JA(J)
              Jb=I_DIAG(JP)
              DO idx=1,lenBlock
                IS=LocalColor % ISindex(iBlock, idx)
                ID=LocalColor % IDindex(iBlock, idx)
                eCoeff=SolDat % ASPAR_pc(IS,ID,J)
!               This is the mystery: the two formulas for eCoeffB gives
!               different results
!                hVal=ONE/SolDat % ASPAR_block(IS,ID,Jb)
!                IF (abs(hVal - ACtest(IS,ID,JP)) .gt. THR8) THEN
!                  WRITE(740+myrank,*) 'hVal=', hVal, 'ACtest=', ACtest(IS,ID,JP)
!                  WRITE(740+myrank,*) '      diff=', hVal - ACtest(IS,ID,JP)
!                  FLUSH(740+myrank)
!                END IF
!                eCoeffB=SolDat % ASPAR_block(IS,ID,J)*hVal
!                eCoeffB=SolDat % ASPAR_block(IS,ID,J)/SolDat % ASPAR_block(IS,ID,Jb)
!                IF (abs(eCoeff - eCoeffB) .gt. THR) THEN
!                  WRITE(740+myrank,*) '1J=', J, 'eCoeff=', eCoeff, 'eCoeffB=', eCoeffB
!                  WRITE(740+myrank,*) '      diff=', eCoeff - eCoeffB
!                  FLUSH(740+myrank)
!                END IF


                ACret(IS,ID,IP)=ACret(IS,ID,IP) - eCoeff*ACret(IS,ID,JP)
              END DO
            END IF
          END DO
# endif
        ENDIF
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_PARTIAL_SOLVE_U(LocalColor, SolDat, iBlock, ACret)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : MNP, IA, JA, I_DIAG, MSC, MDC, rkind, NP_RES, ONE, THR
      USE datapool, only : myrank, iplg
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(in) :: SolDat
      integer, intent(in) :: iBlock
      real(rkind), intent(inout) :: ACret(MSC,MDC,MNP)
      real(rkind) :: eCoeff
      integer lenBlock, IP, JP, idx, J, IS, ID
      lenBlock=LocalColor % BlockLength(iBlock)
      DO IP=NP_RES,1,-1
        IF (LocalColor % CovLower(IP) == 1) THEN
# if defined REORDER_ASPAR_PC
          DO J=LocalColor % IA_U(IP),LocalColor % IA_U(IP+1)-2
            JP=LocalColor % JA_LU(J)
            DO idx=1,lenBlock
              IS=LocalColor % ISindex(iBlock, idx)
              ID=LocalColor % IDindex(iBlock, idx)
              eCoeff=SolDat % ASPAR_pc(IS,ID,J)
              ACret(IS,ID,IP)=ACret(IS,ID,IP) - eCoeff*ACret(IS,ID,JP)
            END DO
          END DO
          J=LocalColor % IA_U(IP+1)-1
          DO idx=1,lenBlock
            IS=LocalColor % ISindex(iBlock, idx)
            ID=LocalColor % IDindex(iBlock, idx)
            ACret(IS,ID,IP)=ACret(IS,ID,IP)*SolDat % ASPAR_pc(IS,ID,J)
          END DO
# elif defined SOR_DIRECT
          DO J=IA(IP),IA(IP+1)-1
            IF (LocalColor % Jstatus_U(J) .eq. 1) THEN
              JP=JA(J)
              DO idx=1,lenBlock
                IS=LocalColor % ISindex(iBlock, idx)
                ID=LocalColor % IDindex(iBlock, idx)
                eCoeff=SolDat % ASPAR_block(IS,ID,J)
                ACret(IS,ID,IP)=ACret(IS,ID,IP) - eCoeff*ACret(IS,ID,JP)
              END DO
            END IF
          END DO
          J=I_DIAG(IP)
          DO idx=1,lenBlock
            IS=LocalColor % ISindex(iBlock, idx)
            ID=LocalColor % IDindex(iBlock, idx)
            ACret(IS,ID,IP)=ACret(IS,ID,IP)/SolDat % ASPAR_block(IS,ID,J)
          END DO
# else
          DO J=IA(IP),IA(IP+1)-1
            IF (LocalColor % Jstatus_U(J) .eq. 1) THEN
              JP=JA(J)
              DO idx=1,lenBlock
                IS=LocalColor % ISindex(iBlock, idx)
                ID=LocalColor % IDindex(iBlock, idx)
                eCoeff=SolDat % ASPAR_pc(IS,ID,J)
!                eCoeffB=SolDat % ASPAR_block(IS,ID,J)
!                IF (abs(eCoeff - eCoeffB) .gt. THR) THEN
!                  WRITE(740+myrank,*) '2J=', J, 'eCoeff=', eCoeff, 'eCoeffB=', eCoeffB
!                  WRITE(740+myrank,*) '      diff=', eCoeff - eCoeffB
!                  FLUSH(740+myrank)
!                END IF
                ACret(IS,ID,IP)=ACret(IS,ID,IP) - eCoeff*ACret(IS,ID,JP)
              END DO
            END IF
          END DO
          J=I_DIAG(IP)
          DO idx=1,lenBlock
            IS=LocalColor % ISindex(iBlock, idx)
            ID=LocalColor % IDindex(iBlock, idx)
!            eCoeff=SolDat % ASPAR_pc(IS,ID,J)
!            eCoeffB=ONE/SolDat % ASPAR_block(IS,ID,J)
!            IF (abs(eCoeff - eCoeffB) .gt. THR) THEN
!              WRITE(740+myrank,*) '3J=', J, 'eCoeff=', eCoeff, 'eCoeffB=', eCoeffB
!              WRITE(740+myrank,*) '      diff=', eCoeff - eCoeffB
!              FLUSH(740+myrank)
!            END IF
            ACret(IS,ID,IP)=ACret(IS,ID,IP)*SolDat % ASPAR_pc(IS,ID,J)
          END DO
# endif
        ENDIF
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_SYNC_SENDRECV(LocalColor, AC)
      USE DATAPOOL, only : LocalColorInfo, MNP, MSC, MDC, rkind
      USE datapool, only : comm, ierr, myrank
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      real(rkind), intent(inout) :: AC(MSC, MDC, MNP)
      INTEGER :: iRank
      integer iSync
      integer nbSync_send, nbSync_recv
      nbSync_recv=LocalColor%sync_nnbr_recv
      nbSync_send=LocalColor%sync_nnbr_send
      DO iSync=1,nbSync_send
        iRank=LocalColor % sync_ListNeigh_send(iSync)
        CALL mpi_isend(AC, 1, LocalColor%sync_p2dsend_type(iSync), iRank, 1009, comm, LocalColor%sync_p2dsend_rqst(iSync), ierr)
      END DO
      DO iSync=1,nbSync_recv
        iRank=LocalColor % sync_ListNeigh_recv(iSync)
        call mpi_irecv(AC,1,LocalColor%sync_p2drecv_type(iSync),iRank,1009,comm,LocalColor%sync_p2drecv_rqst(iSync),ierr)
      END DO
      IF (nbSync_send > 0) THEN
        call mpi_waitall(nbSync_send, LocalColor%sync_p2dsend_rqst, LocalColor%sync_p2dsend_stat,ierr)
      END IF
      IF (nbSync_recv > 0) THEN
        call mpi_waitall(nbSync_recv, LocalColor%sync_p2drecv_rqst, LocalColor%sync_p2drecv_stat,ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_APPLY_PRECOND(LocalColor, SolDat, ACret)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : MSC, MDC, MNP, rkind
      USE datapool, only : myrank
      USE DATAPOOL, only : DO_SYNC_UPP_2_LOW, DO_SYNC_LOW_2_UPP, DO_SYNC_FINAL
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      type(I5_SolutionData), intent(in) :: SolDat
      real(rkind), intent(inout) :: ACret(MSC, MDC, MNP)
      integer iBlock
      DO iBlock=1,LocalColor % Nblock
        IF (DO_SYNC_LOW_2_UPP) THEN
          CALL I5B_EXCHANGE_P3_LOW_2_UPP_Recv(LocalColor, ACret, iBlock)
        END IF
        CALL I5B_PARTIAL_SOLVE_L(LocalColor, SolDat, iBlock, ACret)
        IF (DO_SYNC_LOW_2_UPP) THEN
          CALL I5B_EXCHANGE_P3_LOW_2_UPP_Send(LocalColor, ACret, iBlock)
        END IF
      END DO
      DO iBlock=1,LocalColor%Nblock
        IF (DO_SYNC_UPP_2_LOW) THEN
          CALL I5B_EXCHANGE_P3_UPP_2_LOW_Recv(LocalColor, ACret, iBlock)
        END IF
        CALL I5B_PARTIAL_SOLVE_U(LocalColor, SolDat, iBlock, ACret)
        IF (DO_SYNC_UPP_2_LOW) THEN
          CALL I5B_EXCHANGE_P3_UPP_2_LOW_Send(LocalColor, ACret, iBlock)
        END IF
      END DO
      IF (DO_SYNC_FINAL) THEN
        CALL I5B_SYNC_SENDRECV(LocalColor, ACret)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_APPLY_FCT(LocalColor, SolDat,  ACin, ACret)
      USE DATAPOOL, only : I5_SolutionData, IA, JA, NP_RES, MDC, MSC, MNP, rkind
      USE DATAPOOL, only : LocalColorInfo
      USE datapool, only : exchange_p4d_wwm, myrank
      implicit none
      integer IP, J, idx
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      REAL(rkind), intent(in) :: ACin(MSC, MDC, MNP)
      REAL(rkind), intent(inout) :: ACret(MSC, MDC, MNP)
      REAL(rkind) :: eSum(MSC,MDC)
#ifdef DEBUG
      REAL(rkind) :: Lerror
      REAL(rkind) :: ACtest1(MSC, MDC, MNP)
      REAL(rkind) :: ACtest2(MSC, MDC, MNP)
#endif
      DO IP=1,NP_RES
        eSum=0
        DO J=IA(IP),IA(IP+1)-1
          idx=JA(J)
          eSum=eSum + SolDat % ASPAR_block(:,:,J)*ACin(:,:,idx)
        END DO
        ACret(:,:,IP)=eSum
      END DO
#ifdef DEBUG
      ACtest1=ACret
      ACtest2=ACret
      CALL I5B_EXCHANGE_P4D_WWM(LocalColor, ACtest1)
      CALL EXCHANGE_P4D_WWM(ACtest2)
      !
      CALL I5B_TOTAL_COHERENCY_ERROR_NPRES(ACret, Lerror)
      WRITE(740+myrank,*) 'ACret   NP_RES cohenrency error=', Lerror
      CALL I5B_TOTAL_COHERENCY_ERROR_NPRES(ACtest1, Lerror)
      WRITE(740+myrank,*) 'ACtest1 NP_RES cohenrency error=', Lerror
      CALL I5B_TOTAL_COHERENCY_ERROR_NPRES(ACtest2, Lerror)
      WRITE(740+myrank,*) 'ACtest2 NP_RES cohenrency error=', Lerror
      !
      CALL I5B_TOTAL_COHERENCY_ERROR(ACret, Lerror)
      WRITE(740+myrank,*) 'ACret   MNP coherency error=', Lerror
      CALL I5B_TOTAL_COHERENCY_ERROR(ACtest1, Lerror)
      WRITE(740+myrank,*) 'ACtest1 MNP coherency error=', Lerror
      CALL I5B_TOTAL_COHERENCY_ERROR(ACtest2, Lerror)
      WRITE(740+myrank,*) 'ACtest2 MNP coherency error=', Lerror
      !
      FLUSH(740+myrank)


#endif
# ifdef NO_SELFE_EXCH
      CALL I5B_EXCHANGE_P4D_WWM(LocalColor, ACret)
# else
      CALL EXCHANGE_P4D_WWM(ACret)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE REPLACE_NAN_ZERO(LocalColor, LScal)
      USE DATAPOOL, only : rkind, MSC, MDC, LocalColorInfo
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      real(rkind), intent(inout) :: LScal(MSC,MDC)
      integer IS, ID
      DO IS=1,MSC
        DO ID=1,MDC
          IF (LScal(IS,ID) .ne. LScal(IS,ID)) THEN
            LScal(IS,ID)=0
          END IF
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_L2_LINF(ACw1, ACw2, Norm_L2, Norm_LINF)
      USE DATAPOOL, only : rkind, MNP, MDC, MSC
      USE DATAPOOL, only : nwild_loc_res, NP_RES
      USE datapool, only : myrank, comm, ierr, nproc, istatus, rtype
      implicit none
      real(rkind), intent(in) :: ACw1(MSC, MDC, MNP)
      real(rkind), intent(in) :: ACw2(MSC, MDC, MNP)
      real(rkind), intent(inout) :: Norm_L2(MSC, MDC)
      real(rkind), intent(inout) :: Norm_LINF(MSC, MDC)
      real(rkind) :: LScal(MSC, MDC, 2)
      real(rkind) :: RScal(MSC, MDC, 2)
      integer IP, iProc, IS, ID
      LScal=0
      DO IP=1,NP_RES
        LScal(:,:,1)=LScal(:,:,1) + nwild_loc_res(IP)*((ACw1(:,:,IP) - ACw2(:,:,IP))**2)
        DO IS=1,MSC
          DO ID=1,MDC
            LScal(IS,ID,2)=max(LScal(IS,ID,2), abs(ACw1(IS,ID,IP) - ACw2(IS,ID,IP)))
          END DO
        END DO
      END DO
      IF (myrank == 0) THEN
        DO iProc=2,nproc
          CALL MPI_RECV(RScal,MSC*MDC*2,rtype, iProc-1, 19, comm, istatus, ierr)
          LScal(:,:,1) = LScal(:,:,1) + RScal(:,:,1)
          DO IS=1,MSC
            DO ID=1,MDC
              LScal(IS,ID,2)=max(LScal(IS,ID,2), RScal(IS,ID,2))
            END DO
          END DO
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(LScal,MSC*MDC*2,rtype, iProc-1, 23, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(LScal,MSC*MDC*2,rtype, 0, 19, comm, ierr)
        CALL MPI_RECV(LScal,MSC*MDC*2,rtype, 0, 23, comm, istatus, ierr)
      END IF
      Norm_L2=LScal(:,:,1)
      Norm_LINF=LScal(:,:,2)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_SCALAR(ACw1, ACw2, LScal)
      USE DATAPOOL, only : rkind, MNP, MDC, MSC
      USE DATAPOOL, only : nwild_loc_res, NP_RES
      USE datapool, only : myrank, comm, ierr, nproc, istatus, rtype
      implicit none
      real(rkind), intent(in) :: ACw1(MSC, MDC, MNP)
      real(rkind), intent(in) :: ACw2(MSC, MDC, MNP)
      real(rkind), intent(inout) :: LScal(MSC, MDC)
      real(rkind) :: RScal(MSC, MDC)
      integer IP, iProc
      LScal=0
      DO IP=1,NP_RES
        LScal=LScal + nwild_loc_res(IP)*ACw1(:,:,IP)*ACw2(:,:,IP)
      END DO
      IF (myrank == 0) THEN
        DO iProc=2,nproc
          CALL MPI_RECV(RScal,MSC*MDC,rtype, iProc-1, 19, comm, istatus, ierr)
          LScal = LScal + RScal
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(LScal,MSC*MDC,rtype, iProc-1, 23, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(LScal,MSC*MDC,rtype, 0, 19, comm, ierr)
        CALL MPI_RECV(LScal,MSC*MDC,rtype, 0, 23, comm, istatus, ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_SUM_MAX(ACw, LSum, LMax)
      USE DATAPOOL, only : rkind, MNP, MDC, MSC
      USE DATAPOOL, only : nwild_loc_res, NP_RES
      USE datapool, only : myrank, comm, ierr, nproc, istatus, rtype
      implicit none
      real(rkind), intent(in) :: ACw(MSC, MDC, MNP)
      real(rkind), intent(inout) :: LSum(MSC, MDC)
      real(rkind), intent(inout) :: LMax(MSC, MDC)
      real(rkind) :: RScal(MSC, MDC)
      integer IP, iProc, IS, ID
      LSum=0
      DO IP=1,NP_RES
        LSum=LSum + nwild_loc_res(IP)*ACw(:,:, IP)
      END DO
      DO IS=1,MSC
        DO ID=1,MDC
          LMax(IS,ID)=maxval(ACw(IS,ID,:))
        END DO
      END DO
      IF (myrank == 0) THEN
        DO iProc=2,nproc
          CALL MPI_RECV(RScal,MSC*MDC,rtype, iProc-1, 53, comm, istatus, ierr)
          LSum = LSum + RScal
        END DO
        DO iProc=2,nproc
          CALL MPI_RECV(RScal,MSC*MDC,rtype, iProc-1, 59, comm, istatus, ierr)
          DO IS=1,MSC
            DO ID=1,MDC
              IF (RScal(IS,ID) .gt. LMax(IS,ID)) THEN
                LMax(IS,ID)=RScal(IS,ID)
              END IF
            END DO
          END DO
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(LSum,MSC*MDC,rtype, iProc-1, 197, comm, ierr)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(LMax,MSC*MDC,rtype, iProc-1, 199, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(LSum,MSC*MDC,rtype, 0, 53, comm, ierr)
        CALL MPI_SEND(LMax,MSC*MDC,rtype, 0, 59, comm, ierr)
        CALL MPI_RECV(LSum,MSC*MDC,rtype, 0, 197, comm, istatus, ierr)
        CALL MPI_RECV(LMax,MSC*MDC,rtype, 0, 199, comm, istatus, ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      FUNCTION I5B_SUMTOT(ACw)
      USE DATAPOOL, only : rkind, MNP, MDC, MSC
      implicit none
      real(rkind), intent(in) :: ACw(MSC, MDC, MNP)
      real(rkind) :: LSum(MSC, MDC)
      real(rkind) :: LMax(MSC, MDC)
      real(rkind) :: I5B_SUMTOT
      CALL I5B_SUM_MAX(ACw, LSum, LMax)
      I5B_SUMTOT=sum(LSum)
      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
! We use the notations of
! http://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
! and we use K1=Id
! In this algorithm, the use of v_{i-1}, v_i can be replace to just "v"
! The same for x, r
! 
      SUBROUTINE I5B_BCGS_REORG_SOLVER(LocalColor, SolDat, nbIter, Norm_L2, Norm_LINF)
      USE DATAPOOL, only : MSC, MDC, MNP, NP_RES, NNZ, AC2, WAE_SOLVERTHR
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData, rkind
      USE DATAPOOL, only : PCmethod, STAT, myrank
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      integer, intent(inout) :: nbIter
      REAL(rkind), intent(inout) :: Norm_L2(MSC,MDC)
      REAL(rkind), intent(inout) :: Norm_LINF(MSC,MDC)
      REAL(rkind) :: Rho(MSC,MDC)
      REAL(rkind) :: Prov(MSC,MDC)
      REAL(rkind) :: Alpha(MSC,MDC)
      REAL(rkind) :: Beta(MSC,MDC)
      REAL(rkind) :: Omega(MSC,MDC)
      REAL(rkind) :: MaxError, CritVal
      integer :: MaxIter = 30
      integer IP
      MaxError=WAE_SOLVERTHR
      CALL I5B_APPLY_FCT(LocalColor, SolDat,  AC2, SolDat % AC3)
      SolDat % AC1=0                               ! y
      SolDat % AC3=SolDat % B_block - SolDat % AC3 ! r residual
      SolDat % AC4=SolDat % AC3                    ! hat{r_0} term
      SolDat % AC5=0                               ! v
      SolDat % AC6=0                               ! p
      SolDat % AC7=0                               ! t
      Rho=1
      Alpha=1
      Omega=1
      nbIter=0
!      WRITE(740+myrank,*) 'Beginning solution'
!      FLUSH(740+myrank)
      DO
        nbIter=nbIter+1
!        WRITE(740+myrank,*) 'nbIter=', nbIter
!        FLUSH(740+myrank)

        ! L1: Rhoi =(\hat{r}_0, r_{i-1}
        CALL I5B_SCALAR(SolDat % AC4, SolDat % AC3, Prov)

        ! L2: Beta=(RhoI/Rho(I-1))  *  (Alpha/Omega(i-1))
        Beta=(Prov/Rho)*(Alpha/Omega)
        CALL REPLACE_NAN_ZERO(LocalColor, Beta)
        Rho=Prov

        ! L3: Pi = r(i-1) + Beta*(p(i-1) -omega(i-1)*v(i-1))
        DO IP=1,MNP
          SolDat%AC6(:,:,IP)=SolDat%AC3(:,:,IP)                        &
     &      + Beta(:,:)*SolDat%AC6(:,:,IP)                            &
     &      - Beta(:,:)*Omega(:,:)*SolDat%AC5(:,:,IP)
        END DO

        ! L4 y=K^(-1) Pi
        SolDat%AC1=SolDat%AC6
        IF (PCmethod .gt. 0) THEN
          CALL I5B_APPLY_PRECOND(LocalColor, SolDat, SolDat%AC1)
        ENDIF

        ! L5 vi=Ay
        CALL I5B_APPLY_FCT(LocalColor, SolDat,  SolDat%AC1, SolDat%AC5)

        ! L6 Alpha=Rho/(hat(r)_0, v_i)
        CALL I5B_SCALAR(SolDat % AC4, SolDat % AC5, Prov)
        Alpha(:,:)=Rho(:,:)/Prov(:,:)
        CALL REPLACE_NAN_ZERO(LocalColor, Alpha)

        ! L6.1 x(i)=x(i-1) + Alpha y
        DO IP=1,MNP
          AC2(:,:,IP)=AC2(:,:,IP)                        &
     &      + Alpha(:,:)*SolDat%AC1(:,:,IP)
        END DO

        ! L7 s=r(i-1) - alpha v(i)
        DO IP=1,MNP
          SolDat%AC3(:,:,IP)=SolDat%AC3(:,:,IP)                        &
     &      - Alpha(:,:)*SolDat%AC5(:,:,IP)
        END DO

        ! L8 z=K^(-1) s
        SolDat%AC1=SolDat%AC3
        IF (PCmethod .gt. 0) THEN
          CALL I5B_APPLY_PRECOND(LocalColor, SolDat, SolDat%AC1)
        END IF

        ! L9 t=Az
        CALL I5B_APPLY_FCT(LocalColor, SolDat,  SolDat%AC1, SolDat%AC7)

        ! L10 omega=(t,s)/(t,t)
        CALL I5B_SCALAR(SolDat % AC7, SolDat % AC3, Omega)
        CALL I5B_SCALAR(SolDat % AC7, SolDat % AC7, Prov)
        Omega(:,:)=Omega(:,:)/Prov(:,:)
        CALL REPLACE_NAN_ZERO(LocalColor, Omega)

        ! L11 x(i)=x(i-1) + Omega z
        DO IP=1,MNP
          AC2(:,:,IP)=AC2(:,:,IP)                        &
     &      + Omega(:,:)*SolDat%AC1(:,:,IP)
        END DO

        ! L12 If x is accurate enough finish
        CALL I5B_APPLY_FCT(LocalColor, SolDat,  AC2, SolDat%AC1)
        CALL I5B_L2_LINF(SolDat%AC1, SolDat%B_block, Norm_L2, Norm_LINF)
        CritVal=maxval(Norm_L2)
!        WRITE(740+myrank,*) 'CritVal=', CritVal
!        FLUSH(740+myrank)
        IF (CritVal .lt. MaxError) THEN
          EXIT
        ENDIF
        IF (nbIter .gt. MaxIter) THEN
          EXIT
        ENDIF

        ! L13 r=s-omega t
        DO IP=1,MNP
          SolDat%AC3(:,:,IP)=SolDat%AC3(:,:,IP)                        &
     &      - Omega(:,:)*SolDat%AC7(:,:,IP)
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_ALLOCATE(SolDat)
      USE DATAPOOL, only : I5_SolutionData, MNP, MSC, MDC, NNZ
      implicit none
      type(I5_SolutionData), intent(inout) :: SolDat
      integer istat
      allocate(SolDat % AC1(MSC,MDC,MNP), SolDat % AC3(MSC,MDC,MNP), SolDat % AC4(MSC,MDC,MNP), SolDat % AC5(MSC,MDC,MNP), SolDat % AC6(MSC,MDC,MNP), SolDat % AC7(MSC,MDC,MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 75')
      allocate(SolDat % ASPAR_block(MSC,MDC,NNZ), SolDat % B_block(MSC,MDC, MNP), stat=istat)
# ifndef SOR_DIRECT
      allocate(SolDat % ASPAR_pc(MSC,MDC,NNZ), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 77')
# endif

      IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 78')
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WWM_SOLVER_INIT
      implicit none
      CALL I5B_SOLVER_INIT
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_SOLVER_INIT
      USE DATAPOOL
      implicit none
      NblockFreqDir = NB_BLOCK
      CALL SYMM_INIT_COLORING(MainLocalColor, NblockFreqDir)
# ifdef DEBUG
      WRITE(myrank+740,*) 'After SYMM_INIT_COLORING'
      FLUSH(myrank+740)
# endif
      CALL I5B_ALLOCATE(SolDat)
# ifdef DEBUG
      WRITE(myrank+740,*) 'After I5B_ALLOCATE'
      FLUSH(myrank+740)
# endif
      IF (PCmethod .eq. 2) THEN
!        CALL CREATE_ASPAR_EXCHANGE_ARRAY(LocalColor)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_FREE(SolDat)
      USE DATAPOOL, only : I5_SolutionData
      implicit none
      type(I5_SolutionData), intent(inout) :: SolDat
      deallocate(SolDat % AC1, SolDat % AC3, SolDat % AC4, SolDat % AC5, SolDat % AC6, SolDat % AC7)
      deallocate(SolDat % ASPAR_block, SolDat % B_block, SolDat % ASPAR_pc)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_SUM(AC, eSum)
      USE DATAPOOL, only : MNP, MSC, MDC, rkind, ZERO
      implicit none
      real(rkind), intent(in) :: AC(MNP,MSC,MDC)
      real(rkind), intent(out) :: eSum
      integer IP,IS,ID
      eSum=ZERO
      DO IP=1,MNP
        DO IS=1,MSC
          DO ID=1,MDC
            eSum=eSum + AC(IP,IS,ID)
          END DO
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WWM_SOLVER_EIMPS(LocalColor, SolDat)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData, stat
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'ENTERING WWM_SOLVER_EIMPS'
      FLUSH(STAT%FHNDL)
      CALL I5B_EIMPS(LocalColor, SolDat)
      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'FINISHING WWM_SOLVER_EIMPS'
      FLUSH(STAT%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*
!**********************************************************************
      SUBROUTINE EIMPS_B_BLOCK(U,B)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), intent(out) :: B(MSC, MDC, MNP)
      REAL(rkind), intent(in)  :: U(MSC, MDC, MNP)

      INTEGER :: IP, ID, IS
      INTEGER :: IPGL1, IPREL
 
      REAL(rkind) :: TRIA03
!
# ifdef DEBUG
      WRITE(740+myrank,*) 'Begin of EIMPS_B_BLOCK'
# endif

      DO IP=1,MNP
        DO ID=1,MDC
          B(:,ID,IP) = U(:,ID,IP) * IOBPD(ID,IP)*IOBWB(IP)*IOBDP(IP)*SI(IP)
        ENDDO
      END DO
      IF (LBCWA .OR. LBCSP) THEN
        DO IP = 1, IWBMNP
          IF (LINHOM) THEN
            IPrel=IP
          ELSE
            IPrel=1
          ENDIF
          IPGL1 = IWBNDLC(IP)
          B(:,:,IPGL1) = WBAC(:,:,IPrel)  * SI(IPGL1) ! Overwrite ... 
        END DO
      END IF
# if defined DEBUG
      WRITE(3000+myrank,*)  'sum(B     )=', sum(B)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EIMPS_ASPAR_B_BLOCK_SOURCES_TOTAL(U, ASPAR, B)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), intent(inout) :: ASPAR(MSC, MDC, NNZ)
      REAL(rkind), intent(inout) :: B(MSC, MDC, MNP)
      REAL(rkind), intent(in)    :: U(MSC, MDC, MNP)
      INTEGER :: POS_TRICK(3,2)
      REAL(rkind) :: FL11(MSC,MDC), FL12(MSC,MDC), FL21(MSC,MDC), FL22(MSC,MDC), FL31(MSC,MDC), FL32(MSC,MDC)
      REAL(rkind):: CRFS(MSC,MDC,3), K1(MSC,MDC), KM(MSC,MDC,3), K(MSC,MDC,3), TRIA03
      REAL(rkind):: GTEMP2, DELFL, USFM, FLHAB, LIMFAC
# ifndef NO_MEMORY_CX_CY
      REAL(rkind) :: CX(MSC,MDC,MNP), CY(MSC,MDC,MNP)
# else
      REAL(rkind) :: CXY(2,MSC,MDC,3)
      REAL(rkind)      :: DIFRU, USOC, WVC
# endif
# ifndef SINGLE_LOOP_AMATRIX
      REAL(rkind) :: DELTAL(MSC,MDC,3,MNE)
      REAL(rkind) :: KP(MSC,MDC,3,MNE), NM(MSC,MDC,MNE)
      INTEGER     :: POS
# else
      REAL(rkind) :: DELTAL(MSC,MDC,3)
      REAL(rkind) :: KP(MSC,MDC,3), NM(MSC,MDC)
# endif
      INTEGER :: I1, I2, I3
      INTEGER :: IP, ID, IS, IE
      INTEGER :: I, IPGL1, IPrel
      REAL(rkind) :: DTK(MSC,MDC), TMP3(MSC,MDC)
      REAL(rkind) :: LAMBDA(2,MSC,MDC)
# ifdef DEBUG
      WRITE(740+myrank,*) 'Begin of EIMPS_ASPAR_B_BLOCK'
# endif
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2

# ifndef NO_MEMORY_CX_CY
      CALL CADVXY_VECTOR(CX, CY)
# endif
!
!        Calculate countour integral quantities ...
!
# ifdef DEBUG
      WRITE(740+myrank,*) ' Before MNE loop'
# endif
      ASPAR = 0.0_rkind ! Mass matrix ...
      B     = 0.0_rkind ! Right hand side ...
      DO IE = 1, MNE
# ifndef NO_MEMORY_CX_CY
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        LAMBDA(1,:,:) = ONESIXTH * (CX(:,:,I1) + CX(:,:,I2) + CX(:,:,I3))
        LAMBDA(2,:,:) = ONESIXTH * (CY(:,:,I1) + CY(:,:,I2) + CY(:,:,I3))
        K(:,:,1)  = LAMBDA(1,:,:) * IEN(1,IE) + LAMBDA(2,:,:) * IEN(2,IE)
        K(:,:,2)  = LAMBDA(1,:,:) * IEN(3,IE) + LAMBDA(2,:,:) * IEN(4,IE)
        K(:,:,3)  = LAMBDA(1,:,:) * IEN(5,IE) + LAMBDA(2,:,:) * IEN(6,IE)
        FL11(:,:) = CX(:,:,I2)*IEN(1,IE)+CY(:,:,I2)*IEN(2,IE)
        FL12(:,:) = CX(:,:,I3)*IEN(1,IE)+CY(:,:,I3)*IEN(2,IE)
        FL21(:,:) = CX(:,:,I3)*IEN(3,IE)+CY(:,:,I3)*IEN(4,IE)
        FL22(:,:) = CX(:,:,I1)*IEN(3,IE)+CY(:,:,I1)*IEN(4,IE)
        FL31(:,:) = CX(:,:,I1)*IEN(5,IE)+CY(:,:,I1)*IEN(6,IE)
        FL32(:,:) = CX(:,:,I2)*IEN(5,IE)+CY(:,:,I2)*IEN(6,IE)
# else
        DO I=1,3
          IP = INE(I,IE)
          DO IS=1,MSC
            DO ID=1,MDC
              IF (LSECU .OR. LSTCU) THEN
                CXY(1,IS,ID,I) = CG(IS,IP)*COSTH(ID)+CURTXY(IP,1)
                CXY(2,IS,ID,I) = CG(IS,IP)*SINTH(ID)+CURTXY(IP,2)
              ELSE
                CXY(1,IS,ID,I) = CG(IS,IP)*COSTH(ID)
                CXY(2,IS,ID,I) = CG(IS,IP)*SINTH(ID)
              END IF
              IF (LSPHE) THEN
                CXY(1,IS,ID,I) = CXY(1,IS,ID,I)*INVSPHTRANS(IP,1)
                CXY(2,IS,ID,I) = CXY(2,IS,ID,I)*INVSPHTRANS(IP,2)
              END IF
              IF (LDIFR) THEN
                CXY(1,IS,ID,I) = CXY(1,IS,ID,I)*DIFRM(IP)
                CXY(2,IS,ID,I) = CXY(2,IS,ID,I)*DIFRM(IP)
                IF (LSECU .OR. LSTCU) THEN
                  IF (IDIFFR .GT. 1) THEN
                    WVC = SPSIG(IS)/WK(IS,IP)
                    USOC = (COSTH(ID)*CURTXY(IP,1) + SINTH(ID)*CURTXY(IP,2))/WVC
                    DIFRU = ONE + USOC * (ONE - DIFRM(IP))
                  ELSE
                    DIFRU = DIFRM(IP)
                  END IF
                  CXY(1,IS,ID,I) = CXY(1,IS,ID,I) + DIFRU*CURTXY(IP,1)
                  CXY(2,IS,ID,I) = CXY(2,IS,ID,I) + DIFRU*CURTXY(IP,2)
                END IF
              END IF
            END DO
          END DO
        END DO

        LAMBDA(:,:,:) = ONESIXTH * (CXY(:,:,:,1) + CXY(:,:,:,2) + CXY(:,:,:,3))
        K(:,:,1)  = LAMBDA(1,:,:) * IEN(1,IE) + LAMBDA(2,:,:) * IEN(2,IE)
        K(:,:,2)  = LAMBDA(1,:,:) * IEN(3,IE) + LAMBDA(2,:,:) * IEN(4,IE)
        K(:,:,3)  = LAMBDA(1,:,:) * IEN(5,IE) + LAMBDA(2,:,:) * IEN(6,IE)
        FL11(:,:) = CXY(1,:,:,2)*IEN(1,IE)+CXY(2,:,:,2)*IEN(2,IE)
        FL12(:,:) = CXY(1,:,:,3)*IEN(1,IE)+CXY(2,:,:,3)*IEN(2,IE)
        FL21(:,:) = CXY(1,:,:,3)*IEN(3,IE)+CXY(2,:,:,3)*IEN(4,IE)
        FL22(:,:) = CXY(1,:,:,1)*IEN(3,IE)+CXY(2,:,:,1)*IEN(4,IE)
        FL31(:,:) = CXY(1,:,:,1)*IEN(5,IE)+CXY(2,:,:,1)*IEN(6,IE)
        FL32(:,:) = CXY(1,:,:,2)*IEN(5,IE)+CXY(2,:,:,2)*IEN(6,IE)
# endif
        CRFS(:,:,1) = - ONESIXTH *  (TWO *FL31(:,:) + FL32(:,:) + FL21(:,:) + TWO * FL22(:,:) )
        CRFS(:,:,2) = - ONESIXTH *  (TWO *FL32(:,:) + TWO * FL11(:,:) + FL12(:,:) + FL31(:,:) )
        CRFS(:,:,3) = - ONESIXTH *  (TWO *FL12(:,:) + TWO * FL21(:,:) + FL22(:,:) + FL11(:,:) )
        KM = MIN(ZERO,K)
# ifndef SINGLE_LOOP_AMATRIX
        KP(:,:,:,IE) = MAX(ZERO,K)
        DELTAL(:,:,:,IE) = CRFS(:,:,:) - KP(:,:,:,IE)
        NM(:,:,IE)=ONE/MIN(-THR,KM(:,:,1) + KM(:,:,2) + KM(:,:,3))
# else
        KP(:,:,:) = MAX(ZERO,K)
        DELTAL(:,:,:) = CRFS(:,:,:)- KP(:,:,:)
        NM(:,:)=ONE/MIN(-THR,KM(:,:,1) + KM(:,:,2) + KM(:,:,3))
        TRIA03 = ONETHIRD * TRIA(IE)
        DO I=1,3
          IP=INE(I,IE)
          I1=JA_IE(I,1,IE)
          I2=JA_IE(I,2,IE)
          I3=JA_IE(I,3,IE)
          K1(:,:) =  KP(:,:,I)
          DO ID=1,MDC
            DTK(:,ID) =  K1(:,ID) * DT4A * IOBWB(IP) * IOBPD(ID,IP) * IOBDP(IP)
          END DO
          TMP3(:,:)  =  DTK(:,:) * NM(:,:)
          ASPAR(:,:,I1) =  TRIA03+DTK(:,:)- TMP3(:,:) * DELTAL(:,:,I             ) + ASPAR(:,:,I1)
          ASPAR(:,:,I2) =                 - TMP3(:,:) * DELTAL(:,:,POS_TRICK(I,1)) + ASPAR(:,:,I2)
          ASPAR(:,:,I3) =                 - TMP3(:,:) * DELTAL(:,:,POS_TRICK(I,2)) + ASPAR(:,:,I3)
          DO ID=1,MDC
            B(:,ID,IP)  =  B(:,ID,IP) + U(:,ID,IP) * TRIA03 * IOBWB(IP) * IOBPD(ID,IP) * IOBDP(IP) 
          END DO
        END DO
# endif
      END DO
# ifndef SINGLE_LOOP_AMATRIX
      J     = 0    ! Counter ...
      DO IP = 1, NP_RES
        IF (IOBWB(IP) .EQ. 1 .AND. DEP(IP) .GT. DMIN) THEN
          DO I = 1, CCON(IP)
            J = J + 1
            IE    =  IE_CELL(J)
            POS   =  POS_CELL(J)
            K1(:,:)    =  KP(:,:,POS,IE) ! Flux Jacobian
            TRIA03 = ONETHIRD * TRIA(IE)
            DO ID=1,MDC
              DTK(:,ID)   =  K1(:,ID) * DT4A * IOBPD(ID,IP)
            END DO
            TMP3(:,:)  =  DTK(:,:) * NM(:,:,IE)
            I1    =  POSI(1,J) ! Position of the recent entry in the ASPAR matrix ... ASPAR is shown in fig. 42, p.122
            I2    =  POSI(2,J)
            I3    =  POSI(3,J)
            ASPAR(:,:,I1) =  TRIA03+DTK(:,:)- TMP3(:,:) * DELTAL(:,:,POS             ,IE) + ASPAR(:,:,I1)  ! Diagonal entry
            ASPAR(:,:,I2) =                 - TMP3(:,:) * DELTAL(:,:,POS_TRICK(POS,1),IE) + ASPAR(:,:,I2)  ! off diagonal entries ...
            ASPAR(:,:,I3) =                 - TMP3(:,:) * DELTAL(:,:,POS_TRICK(POS,2),IE) + ASPAR(:,:,I3)
            DO ID=1,MDC
              B(:,ID,IP)     =  B(:,ID,IP) + IOBPD(ID,IP)*TRIA03 * U(:,ID,IP)
            END DO
          END DO
        ELSE
          DO I = 1, CCON(IP)
            J = J + 1
            IE    =  IE_CELL(J)
            TRIA03 = ONETHIRD * TRIA(IE)
            I1    =  POSI(1,J) ! Position of the recent entry in the ASPAR matrix ... ASPAR is shown in fig. 42, p.122
            ASPAR(:,:,I1) =  TRIA03 + ASPAR(:,:,I1)  ! Diagonal entry
            B(:,:,IP)     =  ZERO
          END DO
        END IF
      END DO
# endif
      IF (LBCWA .OR. LBCSP) THEN
        DO IP = 1, IWBMNP
          IF (LINHOM) THEN
            IPrel=IP
          ELSE
            IPrel=1
          ENDIF
          IPGL1 = IWBNDLC(IP)
          ASPAR(:,:,I_DIAG(IPGL1)) = SI(IPGL1) ! Set boundary on the diagonal
          B(:,:,IPGL1)             = WBAC(:,:,IPrel)  * SI(IPGL1)
        END DO
      END IF

      IF (ICOMP .GE. 2 .AND. SMETHOD .GT. 0) THEN
        DO IP = 1, NP_RES
          ASPAR(:,:,I_DIAG(IP)) = ASPAR(:,:,I_DIAG(IP)) + IMATDAA(:,:,IP) * SI(IP) * DT4A !* IOBWB(IP) * IOBDP(IP) ! Add source term to the diagonal
          B(:,:,IP)             = B(:,:,IP) + IMATRAA(:,:,IP) * SI(IP) * DT4A !* IOBWB(IP) * IOBDP(IP) ! Add source term to the right hand side
        END DO
      ENDIF

!      WRITE(*,*) SUM(IMATRAA), SUM(IMATDAA)

# if defined DEBUG
      WRITE(3000+myrank,*)  'sum(ASPAR )=', sum(ASPAR)
      WRITE(3000+myrank,*)  'sum(B     )=', sum(B)
      DO IS=1,MSC
        WRITE(3000+myrank,*) 'IS, sum(ASPAR)=', IS, sum(ASPAR(IS,:,:))
      END DO
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5B_EIMPS(LocalColor, SolDat)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : rkind, MSC, MDC, AC2, MNP, NNZ
      USE DATAPOOL, only : PCmethod, IOBPD, ZERO, STAT
      USE datapool, only : myrank, exchange_p4d_wwm
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      integer nbIter
      real(rkind) :: Norm_L2(MSC,MDC), Norm_LINF(MSC,MDC)
!      real(rkind) :: Lerror
# ifndef ASPAR_B_COMPUTE_BLOCK
      real(rkind) :: U(MNP), ASPAR(NNZ), B(MNP)
# endif
      integer IS, ID, IP

      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'ENTERING I5B_EIMPS'
      FLUSH(STAT%FHNDL)

# ifdef DEBUG
      WRITE(740+myrank,*) 'Begin I5B_EIMPS'
# endif
      CALL EIMPS_ASPAR_B_BLOCK_SOURCES_TOTAL(AC2, SolDat%ASPAR_block, SolDat%B_block)
# ifdef DEBUG
      WRITE(740+myrank,*) 'After ASPAR init'
# endif
!# ifdef NO_SELFE_EXCH
!      CALL I5B_EXCHANGE_P4D_WWM(LocalColor, SolDat % B_block)
!# else
!      CALL I5B_TOTAL_COHERENCY_ERROR_NPRES(SolDat%B_block, Lerror)
!      Print *, 'B_block NP_RES cohenrency error=', Lerror
      CALL EXCHANGE_P4D_WWM(SolDat % B_block)
!# endif
# ifdef DEBUG
      WRITE(740+myrank,*) 'After EXCHANGE_P4D_WWM'
# endif
      CALL I5B_EXCHANGE_ASPAR(LocalColor, SolDat%ASPAR_block)
# ifdef DEBUG
      WRITE(740+myrank,*) 'After I5B_EXCHANGE_ASPAR'
# endif
      CALL I5B_CREATE_PRECOND(LocalColor, SolDat, PCmethod)
# ifdef DEBUG
      WRITE(740+myrank,*) 'After I5B_CREATE_PRECOND'
# endif
      CALL I5B_BCGS_REORG_SOLVER(LocalColor, SolDat, nbIter, Norm_L2, Norm_LINF)
      DO IP=1,MNP
        DO ID=1,MDC
          AC2(:,ID,IP)=MAX(ZERO, AC2(:,ID,IP))*MyREAL(IOBPD(ID,IP))
        END DO
      END DO
      WRITE(STAT%FHNDL,*) 'nbIter=', nbIter, 'L2/LINF=', maxval(Norm_L2), maxval(Norm_LINF)
# ifdef DEBUG
      IF (myrank == 1) THEN
        Write(myrank+591,*) 'Clearing ENDING'
      END IF
# endif
      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'FINISHING I5B_EIMPS'
      FLUSH(STAT%FHNDL)
      END SUBROUTINE
#endif
