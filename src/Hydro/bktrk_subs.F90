!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.


!===============================================================================
!===============================================================================
! SCHISM BACKTRACKING SUBROUTINES
!
! subroutine init_inter_btrack
! subroutine inter_btrack
! subroutine btrack
! subroutine quicksearch
! subroutine intersect2
!
!===============================================================================
!===============================================================================

subroutine init_inter_btrack
!-------------------------------------------------------------------------------
! Initialize data-types for inter-subdomain backtracking.
!-------------------------------------------------------------------------------
!#ifdef USE_MPIMODULE
!  use mpi
!#endif
  use schism_glbl
  use schism_msgp
  implicit none
!#ifndef USE_MPIMODULE
  include 'mpif.h'
!#endif
  integer :: blockl(2),types(2),nmm
#if MPIVERSION==1
  integer :: displ(2),base
#elif MPIVERSION==2
  integer(kind=MPI_ADDRESS_KIND) :: displ(2),base
#endif
  type(bt_type) :: bttmp
!-------------------------------------------------------------------------------

  ! Dimension of inter-subdomain btrack arrays for sending and receiving
  call mpi_allreduce(nsa,nmm,1,itype,MPI_MAX,comm,ierr)
  mxnbt=s1_mxnbt*nmm*nvrt

  ! First part of bt_type is block of 8 integers
  ! (starting at bttmp%rank)
  blockl(1)=8
#if MPIVERSION==1
! displ(1) is the address
  call mpi_address(bttmp%rank,displ(1),ierr)
#elif MPIVERSION==2
  call mpi_get_address(bttmp%rank,displ(1),ierr)
#endif
  if(ierr/=MPI_SUCCESS) call parallel_abort('INIT_INTER_BTRACK: mpi_get_address',ierr)
  types(1)=itype

  ! Second part of bt_type is block of 26+mntracers doubles (including arrays)
  ! (starting at bttmp%dtbk)
  blockl(2)=26+mntracers
#if MPIVERSION==1
  call mpi_address(bttmp%dtbk,displ(2),ierr)
#elif MPIVERSION==2
  call mpi_get_address(bttmp%dtbk,displ(2),ierr)
#endif
  if(ierr/=MPI_SUCCESS) call parallel_abort('INIT_INTER_BTRACK: mpi_get_address',ierr)
  types(2)=rtype

  ! Shift displ to compute actual displacements
  base=displ(1)
  displ(1)=displ(1)-base
  displ(2)=displ(2)-base

  ! MPI datatype for bt_type
#if MPIVERSION==1
  call mpi_type_struct(2,blockl,displ,types,bt_mpitype,ierr)
#elif MPIVERSION==2
  call mpi_type_create_struct(2,blockl,displ,types,bt_mpitype,ierr)
#endif
  if(ierr/=MPI_SUCCESS) call parallel_abort('INIT_INTER_BTRACK: type_create',ierr)
  call mpi_type_commit(bt_mpitype,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('INIT_INTER_BTRACK: type_commit',ierr)

end subroutine init_inter_btrack

!===============================================================================
!===============================================================================

subroutine inter_btrack(itime,nbt,btlist)
!-------------------------------------------------------------------------------
! Routine for completing inter-subdomain backtracking.
!
! Input:
!   itime: global time stepping # (info only);
!   nbt: number of inter-subdomain trajectories
!   btlist: list of inter-subdomain trajectories
!
! Output:
!   btlist: list of completed inter-subdomain trajectories; not in original order
!-------------------------------------------------------------------------------
!#ifdef USE_MPIMODULE
!  use mpi
!#endif
  use schism_glbl
  use schism_msgp
  implicit none
!#ifndef USE_MPIMODULE
  include 'mpif.h'
!#endif

  integer,intent(in) :: itime
  integer,intent(in) :: nbt
  type(bt_type),intent(inout) :: btlist(mxnbt)

  integer :: stat,i,ii,j,ie,irank,nnbrq,inbr,nbts,nbtd
  integer :: mxbtsend,mxbtrecv,mnbt
  real(rkind) :: xt,yt,zt,uuint,vvint,wwint,ttint,ssint
  logical :: lexit,bt_donel(1),bt_done(1)
  integer :: icw
  real(rkind) :: cwtmp

  integer :: ncmplt,icmplt(nproc)
  integer,allocatable :: nbtsend(:),ibtsend(:,:),nbtrecv(:),ibtrecv(:,:)
#if MPIVERSION==1
  integer,allocatable :: bbtsend(:),bbtrecv(:)
#endif
  integer,allocatable :: btsend_type(:),btsend_rqst(:),btsend_stat(:,:)
  integer,allocatable :: btrecv_type(:),btrecv_rqst(:),btrecv_stat(:,:)
  type(bt_type),allocatable :: btsendq(:),btrecvq(:),bttmp(:),btdone(:)
#ifdef DEBUG
    integer,save :: ncalls=0
#endif
!-------------------------------------------------------------------------------

#ifdef DEBUG
  ncalls=ncalls+1
  fdb='interbtrack_0000'
  lfdb=len_trim(fdb)
  write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
  if(ncalls==1) then
    open(30,file=out_dir(1:len_out_dir)//fdb,status='replace')
  else
    open(30,file=out_dir(1:len_out_dir)//fdb,status='old',position='append')
  endif
  write(30,'(a,3i6)') 'INTER_BTRACK START: ',itime,nbt
#endif

  ! Index of wall-timer
!  if(imode==1) then
    icw=4
!  else
!    icw=9
!  endif

  ! Compute max nbt for dimension parameter
#ifdef INCLUDE_TIMING
  cwtmp=mpi_wtime()
#endif
  call mpi_allreduce(nbt,mnbt,1,itype,MPI_MAX,comm,ierr) !mnbt>=1
#ifdef INCLUDE_TIMING
  wtimer(4,2)=wtimer(4,2)+mpi_wtime()-cwtmp
#endif
  mnbt=mnbt*s2_mxnbt !add a scale
 
  ! Allocate type bt_type arrays for sending and receiving
  allocate(btsendq(mnbt),btrecvq(mnbt*nnbr),bttmp(mnbt),btdone(mnbt),stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: type bt_type allocation failure')

  ! Allocate communication data structures
  allocate(nbtsend(nnbr),ibtsend(mnbt,nnbr), &
  nbtrecv(nnbr),ibtrecv(mnbt,nnbr), &
  btsend_type(nnbr),btsend_rqst(nnbr),btsend_stat(MPI_STATUS_SIZE,nnbr), &
  btrecv_type(nnbr),btrecv_rqst(nnbr),btrecv_stat(MPI_STATUS_SIZE,nnbr),stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: comm data allocation failure')
#if MPIVERSION==1
  allocate(bbtsend(mnbt),bbtrecv(mnbt),stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: bbtsend/recv allocation failure')
  bbtsend=1; bbtrecv=1; !blocksize is always 1
#endif

!  write(30,*)'Mark 1',nnbr,mxnbt,iegrpv(btlist(1:nbt)%iegb)

  ! Initialize send queue
  nbts=nbt
  btsendq(1:nbts)=btlist(1:nbts)

  ! Initialize completed bt count
  nbtd=0

  !-----------------------------------------------------------------------------
  ! Outer loop:
  ! > All ranks participate until all inter-subdomain backtracked trajectories
  !   are completed.
  ! > Completed trajectories are placed in btdone list.
  !-----------------------------------------------------------------------------
  outer_loop: do

#ifdef INCLUDE_TIMING
  ! Init communication timer
  cwtmp=mpi_wtime()
#endif

!  write(30,*)'Mark 2'

  ! Count and index sends
  nbtsend=0
  do i=1,nbts
    irank=iegrpv(btsendq(i)%iegb)
    inbr=ranknbr(irank)
    if(inbr==0) then
      write(errmsg,*) 'INTER_BTRACK: bt to non-neighbor!',irank
      call parallel_abort(errmsg)
    endif
    nbtsend(inbr)=nbtsend(inbr)+1
    if(nbtsend(inbr)>mnbt) call parallel_abort('bktrk_subs: overflow (1)')
    ibtsend(nbtsend(inbr),inbr)=i-1 !displacement
  enddo

  ! Set MPI bt send datatypes
  do inbr=1,nnbr
    if(nbtsend(inbr)/=0) then
#if MPIVERSION==1
      call mpi_type_indexed(nbtsend(inbr),bbtsend,ibtsend(1,inbr),bt_mpitype, &
      btsend_type(inbr),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nbtsend(inbr),1,ibtsend(1,inbr),bt_mpitype, &
      btsend_type(inbr),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: create btsend_type',ierr)
      call mpi_type_commit(btsend_type(inbr),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: commit btsend_type',ierr)
    endif
  enddo !inbr

  ! Post recvs for bt counts
  do inbr=1,nnbr
    call mpi_irecv(nbtrecv(inbr),1,itype,nbrrank(inbr),700,comm,btrecv_rqst(inbr),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: irecv 700',ierr)
  enddo

  ! Post sends for bt counts
  do inbr=1,nnbr
    call mpi_isend(nbtsend(inbr),1,itype,nbrrank(inbr),700,comm,btsend_rqst(inbr),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: isend 700',ierr)
  enddo

!  write(30,*)'Mark 4'

  ! Wait for recvs to complete
  call mpi_waitall(nnbr,btrecv_rqst,btrecv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: waitall recv 700',ierr)
  ! Wait for sends to complete
  call mpi_waitall(nnbr,btsend_rqst,btsend_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: waitall send 700',ierr)

!  write(30,*)'Mark 4.5'

  ! Set MPI bt recv datatypes
  i=0 !total # of recv from all neighbors 
  nnbrq=0; !# of "active" neighbors
  do inbr=1,nnbr
    if(nbtrecv(inbr)/=0) then
      nnbrq=nnbrq+1
      if(nbtrecv(inbr)>mnbt) call parallel_abort('bktrk_subs: overflow (3)')
      do j=1,nbtrecv(inbr); ibtrecv(j,inbr)=i+j-1; enddo; !displacement
      i=i+nbtrecv(inbr)
#if MPIVERSION==1
      call mpi_type_indexed(nbtrecv(inbr),bbtrecv,ibtrecv(1,inbr),bt_mpitype, &
      btrecv_type(inbr),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nbtrecv(inbr),1,ibtrecv(1,inbr),bt_mpitype, &
      btrecv_type(inbr),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: create btrecv_type',ierr)
      call mpi_type_commit(btrecv_type(inbr),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: commit btrecv_type',ierr)
!      write(30,*)'recev from',nbrrank(inbr),nbtrecv(inbr)
    endif
  enddo !inbr
 
  ! Check bound for btrecvq
  if(i>mnbt*nnbr) call parallel_abort('bktrk_subs: overflow (2)')

!  write(30,*)'Mark 5'

  ! Post sends for bt data
  do inbr=1,nnbr
    if(nbtsend(inbr)/=0) then
      ! btsendq(1)%rank is the starting address
      call mpi_isend(btsendq(1)%rank,1,btsend_type(inbr),nbrrank(inbr),701, &
      comm,btsend_rqst(inbr),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: isend 701',ierr)
    else
      btsend_rqst(inbr)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post recvs for bt data
  do inbr=1,nnbr
    if(nbtrecv(inbr)/=0) then
      call mpi_irecv(btrecvq(1)%rank,1,btrecv_type(inbr),nbrrank(inbr),701, &
      comm,btrecv_rqst(inbr),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: irecv 701',ierr)
    else
      btrecv_rqst(inbr)=MPI_REQUEST_NULL
    endif
  enddo

!  write(30,*)'Mark 6'

#ifdef INCLUDE_TIMING
  ! Add to communication timer
  wtimer(icw,2)=wtimer(icw,2)+mpi_wtime()-cwtmp
#endif

  !-----------------------------------------------------------------------------
  ! Inner loop: (for efficiency)
  ! > Process inter-subdomain backtracking receive queue until empty
  ! > Completed trajectories placed in btdone list
  ! > Exited trajectories placed in "temporary" send queue bttmp()
  !-----------------------------------------------------------------------------
  nbts=0 !current # of requests from myrank to all neighbors for inter-domain tracking
  inner_loop: do
  if(nnbrq==0) exit inner_loop

#ifdef DEBUG
  write(30,'(a)') 'INNER LOOP'
#endif

  ! Wait for some bt recvs to complete
  ! Parameters of mpi_waitsome:
  ! Inputs: nnbr - # of requests (dimension of btrecv_rqst); btrecv_rqst - array of requests;
  ! Outputs: ncmplt - # of completed requests; icmplt - array of indices of completed operations (integer);
  !          btrecv_stat - array of status objects for completed operations.

#ifdef INCLUDE_TIMING
  cwtmp=mpi_wtime()
#endif

  call mpi_waitsome(nnbr,btrecv_rqst,ncmplt,icmplt,btrecv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: waitsome recv 701',ierr)
#ifdef INCLUDE_TIMING
  wtimer(icw,2)=wtimer(icw,2)+mpi_wtime()-cwtmp
#endif

!  write(30,*)'Mark 7'

  ! Perform local backtracking for received trajectories
  ! Process several neighbors at a time as soon as the info is received from them
!$OMP parallel if(ncmplt>nthreads) default(shared) private(ii,inbr,j,i,ie,lexit)
!$OMP do
  do ii=1,ncmplt
    inbr=icmplt(ii)
    do j=1,nbtrecv(inbr)
      i=ibtrecv(j,inbr)+1
      ie=iegl(btrecvq(i)%iegb)%id

!      write(12,*)'btrack #',ii,btrecvq(i) !ielg(ie),btrecvq(i)%jvrt,btrecvq(i)%rank

      call btrack(btrecvq(i)%l0,btrecvq(i)%i0gb,btrecvq(i)%isbndy,btrecvq(i)%j0, &
&btrecvq(i)%adv,btrecvq(i)%gcor0,btrecvq(i)%frame0,btrecvq(i)%dtbk,btrecvq(i)%vis, &
&btrecvq(i)%rt,btrecvq(i)%rt2,btrecvq(i)%ut,btrecvq(i)%vt,btrecvq(i)%wt, &
&ie,btrecvq(i)%jvrt,btrecvq(i)%xt,btrecvq(i)%yt,btrecvq(i)%zt, &
&btrecvq(i)%sclr,lexit)

!      write(12,*)'btrack #',ii,btrecvq(i) !ielg(ie),btrecvq(i)%jvrt,btrecvq(i)%rank

      if(lexit) then !backtracking exits augmented subdomain
!$OMP   critical
        !Move point to "temporary" send queue
        nbts=nbts+1
        btrecvq(i)%iegb=ielg(ie)
        if(nbts>mnbt) call parallel_abort('bktrk_subs: overflow (5)')
        bttmp(nbts)=btrecvq(i)
!$OMP   end critical
      else !backtracking completed within augmented subdomain
!$OMP   critical
        !Move point to done list
        nbtd=nbtd+1
        btrecvq(i)%iegb=ielg(ie)
        if(nbtd>mnbt) call parallel_abort('bktrk_subs: overflow (4)')
        btdone(nbtd)=btrecvq(i)
!$OMP   end critical
      endif
    enddo !j
  enddo !ii=1,ncmplt
!$OMP end do
!$OMP end parallel

  ! Decrement nnbrq according to number of recvs processed
  nnbrq=nnbrq-ncmplt

  !-----------------------------------------------------------------------------
  ! End inner loop
  !-----------------------------------------------------------------------------
  enddo inner_loop
 
#ifdef DEBUG
  write(30,'(a)') 'DONE INNER LOOP'
#endif

  ! All recv's are complete

#ifdef INCLUDE_TIMING
  ! Init communication timer
  cwtmp=mpi_wtime()
#endif

  ! Wait for sends to complete
  call mpi_waitall(nnbr,btsend_rqst,btsend_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: waitall send 701',ierr)

  ! Free MPI bt send datatypes
  do inbr=1,nnbr
    if(nbtsend(inbr)/=0) then
      call mpi_type_free(btsend_type(inbr),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: free btsend_type',ierr)
    endif
  enddo !inbr

  ! Free MPI bt recv datatypes
  do inbr=1,nnbr
    if(nbtrecv(inbr)/=0) then
      call mpi_type_free(btrecv_type(inbr),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: free btrecv_type',ierr)
    endif
  enddo !inbr

#ifdef DEBUG
  write(30,'(a,4i6)') 'CYCLE OUTER: ',itime,nbts,nbtd
#endif

  ! Exit outer loop if all backtracks are completed
  bt_donel(1)=(nbts==0)
  call mpi_allreduce(bt_donel,bt_done,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
#ifdef INCLUDE_TIMING
  ! Add to communication timer
  wtimer(icw,2)=wtimer(icw,2)+mpi_wtime()-cwtmp
#endif
  if(bt_done(1)) exit outer_loop

  ! Initialize send queue for next cycle of outer loop
  btsendq(1:nbts)=bttmp(1:nbts)


  !-----------------------------------------------------------------------------
  ! End outer loop
  !-----------------------------------------------------------------------------
  enddo outer_loop

  ! Deallocate some type bt_type arrays
  deallocate(btsendq,btrecvq,bttmp,stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: bt type deallocation failure (1)')

  !-----------------------------------------------------------------------------
  ! All inter-subdomain backtracked trajectories completed. Communicate
  ! completed trajectories back to originating subdomain. Note that originating
  ! subdomain is not necessarily a neighbor to the subdomain where the
  ! trajectory was completed.
  ! Avoid communicating to itself as trajectory can come back!
  !-----------------------------------------------------------------------------

#ifdef INCLUDE_TIMING
  ! Init communication timer
  cwtmp=mpi_wtime()
#endif

#ifdef DEBUG
  write(30,'(a)') 'Start all-rank communication'
#endif

  ! Deallocate communication data structures
  deallocate(nbtsend,ibtsend,nbtrecv,ibtrecv, &
  btsend_type,btsend_rqst,btsend_stat, &
  btrecv_type,btrecv_rqst,btrecv_stat)
#if MPIVERSION==1
  deallocate(bbtsend,bbtrecv)
#endif

  ! nbtsend: # of sends from myrank to each proc excluding myrank
  allocate(nbtsend(0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: nbtsend allocation failure')
  nbtsend=0
  do i=1,nbtd
    irank=btdone(i)%rank !destination rank
    if(irank/=myrank) nbtsend(irank)=nbtsend(irank)+1
  enddo !i
  mxbtsend=0; do irank=0,nproc-1; mxbtsend=max(mxbtsend,nbtsend(irank)); enddo;

  ! ibtsend: displacement or position in btdone(1:nbtd) list
  allocate(ibtsend(mxbtsend,0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: ibtsend allocation failure')
#if MPIVERSION==1
  allocate(bbtsend(mxbtsend),stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: bbtsend allocation failure')
  bbtsend=1 !blocksize is always 1
#endif
  nbtsend=0
  do i=1,nbtd
    irank=btdone(i)%rank
    if(irank/=myrank) then
      nbtsend(irank)=nbtsend(irank)+1
      ibtsend(nbtsend(irank),irank)=i-1 !displacement
    endif
  enddo !i

#ifdef DEBUG
  write(30,'(a,66i6)') 'INTER_BTRACK -- NBTSEND: ', &
  itime,(nbtsend(irank),irank=0,nproc-1)
#endif

  ! Allocate recv count array
  ! nbtrecv: # of recv's to myrank from each proc excluding myrank
  allocate(nbtrecv(0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: nbtrecv allocation failure')

  ! All to all scatter/gather of send counts
  call mpi_alltoall(nbtsend,1,itype,nbtrecv,1,itype,comm,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: alltoall nbtsend',ierr)

#ifdef DEBUG
  write(30,'(a,66i6)') 'INTER_BTRACK -- NBTRECV: ', &
  itime,(nbtrecv(irank),irank=0,nproc-1)
#endif

  ! Count max number of completed recv trajectories per rank and allocate ibtrecv
  mxbtrecv=0; do irank=0,nproc-1; mxbtrecv=max(mxbtrecv,nbtrecv(irank)); enddo;
  allocate(ibtrecv(mxbtrecv,0:nproc-1),stat=stat) !position in blist
  if(stat/=0) call parallel_abort('INTER_BTRACK: ibtrecv allocation failure')
#if MPIVERSION==1
  allocate(bbtrecv(mxbtrecv),stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: bbtrecv allocation failure')
  bbtrecv=1 !blocksize is always 1
#endif

  ! Allocate remaining communication data structures
  allocate(btsend_type(0:nproc-1),btsend_rqst(0:nproc-1), &
  btsend_stat(MPI_STATUS_SIZE,0:nproc-1), &
  btrecv_type(0:nproc-1),btrecv_rqst(0:nproc-1), &
  btrecv_stat(MPI_STATUS_SIZE,0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: bt type/rqst/stat allocation failure')

  ! Set MPI bt send datatypes
  do irank=0,nproc-1
    if(nbtsend(irank)/=0) then
      if(irank==myrank) call parallel_abort('INTER_BTRACK: self communication (1)')
#if MPIVERSION==1
      call mpi_type_indexed(nbtsend(irank),bbtsend,ibtsend(1,irank),bt_mpitype, &
      btsend_type(irank),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nbtsend(irank),1,ibtsend(1,irank),bt_mpitype, &
      btsend_type(irank),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: create btsend_type',ierr)
      call mpi_type_commit(btsend_type(irank),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: commit btsend_type',ierr)
    endif
  enddo !irank

  ! Set MPI bt recv datatypes
  i=0 !total # of recv's
  do irank=0,nproc-1
    if(nbtrecv(irank)/=0) then
      if(irank==myrank) call parallel_abort('INTER_BTRACK: self communication (2)')
      do j=1,nbtrecv(irank); ibtrecv(j,irank)=i+j-1; enddo; !displacement
      i=i+nbtrecv(irank)
#if MPIVERSION==1
      call mpi_type_indexed(nbtrecv(irank),bbtrecv,ibtrecv(1,irank),bt_mpitype, &
      btrecv_type(irank),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nbtrecv(irank),1,ibtrecv(1,irank),bt_mpitype, &
      btrecv_type(irank),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: create btrecv_type',ierr)
      call mpi_type_commit(btrecv_type(irank),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: commit btrecv_type',ierr)
    endif
  enddo !irank

  ! Check bound for btlist
  if(i>mxnbt) call parallel_abort('INTER_BTRACK: overflow (6)')  

  ! Post sends for bt data
  do irank=0,nproc-1
    if(nbtsend(irank)/=0) then
      ! irank is checked to not be myrank
      call mpi_isend(btdone(1)%rank,1,btsend_type(irank),irank,711, &
      comm,btsend_rqst(irank),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: isend 711',ierr)
    else
      btsend_rqst(irank)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post recvs for bt data
  do irank=0,nproc-1
    if(nbtrecv(irank)/=0) then
      ! irank is checked to not be myrank
      call mpi_irecv(btlist(1)%rank,1,btrecv_type(irank),irank,711, &
      comm,btrecv_rqst(irank),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: irecv 711',ierr)
    else
      btrecv_rqst(irank)=MPI_REQUEST_NULL
    endif
  enddo

  ! Wait for all bt recvs to complete
  call mpi_waitall(nproc,btrecv_rqst,btrecv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: waitall recv 711',ierr)

  ! Wait for all bt sends to complete
  call mpi_waitall(nproc,btsend_rqst,btsend_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: waitall send 711',ierr)

  ! Append btlist with left-over from btdone list
  do ii=1,nbtd
    irank=btdone(ii)%rank
    if(irank==myrank) then
      i=i+1
      if(i>mxnbt) call parallel_abort('INTER_BTRACK: overflow (7)')  
      btlist(i)=btdone(ii)
#ifdef DEBUG
      write(30,*)'Back to myself!'
#endif
    endif
  enddo !ii

  ! Check if the total # of recv's is nbt
  if(i/=nbt) call parallel_abort('bktrk_subs: mismatch (1)')

  ! Free MPI bt send datatypes
  do irank=0,nproc-1
    if(nbtsend(irank)/=0) then
      call mpi_type_free(btsend_type(irank),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: free btsend_type',ierr)
    endif
  enddo !irank

  ! Free MPI bt recv datatypes
  do irank=0,nproc-1
    if(nbtrecv(irank)/=0) then
      call mpi_type_free(btrecv_type(irank),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: free btrecv_type',ierr)
    endif
  enddo !irank

  ! Deallocate btdone
  deallocate(btdone,stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: btdone deallocation failure')

  ! Deallocate communication data structures
  deallocate(nbtsend,ibtsend,nbtrecv,ibtrecv, &
  &btsend_type,btsend_rqst,btsend_stat, &
  &btrecv_type,btrecv_rqst,btrecv_stat,stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: rqst/stat deallocation failure')
#if MPIVERSION==1
  deallocate(bbtsend,bbtrecv,stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: bbtsend/recv deallocation failure')
#endif

#ifdef INCLUDE_TIMING
  ! Add to communication timer
  wtimer(icw,2)=wtimer(icw,2)+mpi_wtime()-cwtmp
#endif

#ifdef DEBUG
  close(99)
#endif

end subroutine inter_btrack

!===============================================================================
!===============================================================================

      subroutine btrack(l_ns,ipsgb,ifl_bnd,j0,iadvf,gcor0,frame0,dtbk, &
&vis_coe,time_rm,time_rm2,uuint,vvint,wwint,nnel,jlev,xt,yt,zt,sclr,iexit)
!
!***************************************************************************
!									   
! Routine for backtracking. 			   
! For ics=2, tracking is done in the starting frame/plane. 
!       Input:
!             l_ns: side if l_ns=2; element if l_ns=3;
!             ipsgb: global originating node or side or element index (info only);
!             ifl_bnd: Flag for originating node or side on the boundary (for Kriging);
!             j0: Originating vertical level (info only);
!             iadvf: advection flag (0 or 1: use Euler tracking; 2: R-K tracking) 
!                    associated with originating position;
!             gcor0(3): global coord. of the originating pt (for ics=2);
!             frame0(3,3): frame tensor at original start pt (node/side/element) 
!                          (2nd index is axis id) for ics=2;
!             dtbk: target tracking step (actual step may be smaller)
!             vis_coe: weighting factor between continuous and discontinuous vel. - not used anymore
!             
!       Input/output:
!             time_rm: remaining time (<=dt);
!             time_rm2: remaining time (<=dtbk) - leftover from previous subdomain;
!             uuint,vvint,wwint: starting and interpolated vel.; if ics=2, it's
!                                in frame0
!             xt,yt,zt: starting and current coordinates (in local frame0 if ics=2);
!             nnel,jlev: initial and final element and level;
!
!       Output:
!             sclr(4+mntracers): btrack'ed values of some variables;
!             iexit: logical flag indicating backtracking exits augmented subdomain. If
!                    iexit=.true., nnel is inside the aug. domain and should also be inside
!                    one of the neighboring processes. (xt,yt) is inside nnel.
!									   
!***************************************************************************
!
      use schism_glbl
      use schism_msgp, only : parallel_abort,myrank
      implicit none
!      real(rkind), parameter :: per1=1.e-3 !used to check error in tracking
!      real(rkind), parameter :: safety=0.8 !safyty factor in estimating time step

      integer, intent(in) :: l_ns,ipsgb,ifl_bnd,j0,iadvf
      real(rkind), intent(in) :: gcor0(3),frame0(3,3),dtbk,vis_coe
      real(rkind), intent(inout) :: time_rm,time_rm2,uuint,vvint,wwint,xt,yt,zt
      integer, intent(inout) :: nnel,jlev
      real(rkind), intent(out) :: sclr(4+mntracers)
      logical, intent(out) :: iexit

      !Local
      !Function
      real(rkind) :: covar,signa1

      integer :: idt,iflqs1,kbpl,iadptive,nnel0,jlev0,ie,npp,nd,ifl,i,n1, &
                 &n2,n3,kbb,ibelow,isd,in1,in2,j,jj
      real(rkind) :: x0,y0,z0,dtb,trm,zrat,uuint1,vvint1,wwint1,vmag,rr, &
                     &covar2,xn1,xn2,xn3,yn1,yn2,yn3,tmp, &
                     &aa1,aa2,aa3,tnd_min,tnd_max,snd_min,snd_max, &
                     &s_min,s_max,dfvint,dfhint

      real(rkind) :: vxl(3,2),vyl(3,2),vzl(3,2) !,vxn(3),vyn(3),vzn(3)
      real(rkind) :: arco(4),t_xi(6),s_xi(6),sig(3),subrat(4),ztmp(nvrt), &
     &swild(10),swild2(10,nvrt),swild3(nvrt) 
      real(rkind) :: al_beta(mnei_kr+3,4),uvdata(mnei_kr,3) 
      logical :: lrk

!     Constants used in 5th order R-K
!      a(2)=0.2; a(3)=0.3; a(4)=0.6; a(5)=1; a(6)=0.875; a(7)=1
!      b(2,1)=0.2; b(3,1)=0.075; b(3,2)=0.225; b(4,1)=0.3; b(4,2)=-0.9; b(4,3)=1.2
!      b(5,1)=-11./54; b(5,2)=2.5; b(5,3)=-70./27; b(5,4)=35./27
!      b(6,1)=1631./55296; b(6,2)=175./512; b(6,3)=575./13824; b(6,4)=44275./110592; b(6,5)=253./4096
!!     b(7,*) are c(*)
!      b(7,1)=37./378; b(7,2)=0; b(7,3)=250./621; b(7,4)=125./594; b(7,5)=0; b(7,6)=512./1771
!!     dc() are c_i-c_i^*
!      dc(1)=b(7,1)-2825./27648; dc(2)=0; dc(3)=b(7,3)-18575./48384
!      dc(4)=b(7,4)-13525./55296; dc(5)=b(7,5)-277./14336; dc(6)=b(7,6)-0.25

      sclr=0 
      !Initial exit is false
      iexit=.false.

!...  Euler tracking (including left-over from RK2; in this case do Euler only once)
      lrk=iadvf==2.and.time_rm2>0._rkind

      if(iadvf==1.or.iadvf==0.or.lrk) then
!------------------------------------------------------------------------------------------------
      x0=xt
      y0=yt
      z0=zt
      idt=0 !iteration #
      do 
        idt=idt+1
        if(time_rm2>0) then
          !Finish left-over
          dtb=time_rm2
          time_rm2=-99._rkind !reset flag
        else !normal
          dtb=min(dtbk,time_rm)
        endif
        if(dtb<=0._rkind) then
          write(errmsg,*)'BTRACK: dtb<=0,',dtb,dtbk,time_rm,time_rm2,lrk,iadvf
          call parallel_abort(errmsg)
        endif
        xt=x0-dtb*uuint
        yt=y0-dtb*vvint
        zt=z0-dtb*wwint

        call quicksearch(1,idt,l_ns,ipsgb,gcor0,frame0,dtb,x0,y0,z0,nnel,jlev, &
     &xt,yt,zt,trm,iflqs1,kbpl,arco,zrat,ztmp,vis_coe,uuint,vvint,wwint, &
     &uuint1,vvint1,wwint1)

!       Check aug. exit
        if(iflqs1==2) then !exit upon entry into quicksearch
          if(time_rm2>0._rkind) call parallel_abort('BTRACK: just in')
!         Pt (xt,yt,zt) reset, which is inside nnel. 
!         jlev, uuint,vvint,wwint are unchanged (to avoid CPU dependency)
          time_rm2=-99._rkind
          iexit=.true.; return
        endif !iflqs1==2

        if(iflqs1==3) then 
!         Exit during iteration in quicksearch; make sure vel. pointing away from nnel
          if(trm<=0._rkind) call parallel_abort('BTRACK: trm<=0')
          time_rm2=trm
          time_rm=time_rm-(dtb-trm)
          iexit=.true.; return
        endif

        uuint=uuint1; vvint=vvint1; wwint=wwint1
!       Check if vel. is too small to continue
!       This avoids 0 vel case (e.g. hits the no-slip bottom)
        vmag=sqrt(uuint*uuint+vvint*vvint)
        if(vmag<=velmin_btrack.or.iflqs1==1) then
          time_rm=0._rkind; exit
        endif

        if(lrk) then !left-over from RK2; only once
          time_rm=time_rm-dtb
          exit
        endif !lrk

!       Update time_rm
        time_rm=time_rm-(dtb-trm)
        if(time_rm<=real(1.e-6,rkind)*dt) exit

        x0=xt
        y0=yt
        z0=zt
      end do !idt
!------------------------------------------------------------------------------------------------
      endif !Euler

!...  2nd-order R-K tracking
!     Results may vary with # of processors due to crossing of aug. domain
!     but at least error ->0 as more CPUs are used

      if(iadvf==2.and.time_rm>=real(1.e-6,rkind)*dt) then
!------------------------------------------------------------------------------------------------
!     (xt,yt,zt),(uuint etc),nnel,jlev may be carried over from Euler
      x0=xt
      y0=yt
      z0=zt
      iadptive=0 !# of times when dtb is reduced
      idt=0 !iteration #
      dtb=min(dtbk,time_rm) !init.
      if(dtb<=0._rkind) call parallel_abort('BTRACK: dtb<=0 (2a)')
      nnel0=nnel; jlev0=jlev !save 
      do 
        idt=idt+1
        xt=x0-0.5_rkind*dtb*uuint
        yt=y0-0.5_rkind*dtb*vvint
        zt=z0-0.5_rkind*dtb*wwint
        nnel=nnel0; jlev=jlev0 !reset
        call quicksearch(2,idt,l_ns,ipsgb,gcor0,frame0,0.5_rkind*dtb,x0,y0,z0,nnel,jlev, &
     &xt,yt,zt,trm,iflqs1,kbpl,arco,zrat,ztmp,vis_coe,uuint,vvint,wwint, &
     &uuint1,vvint1,wwint1)

        if(iflqs1==2) then !exit upon entry into quicksearch
!         Pt (xt,yt,zt) reset, which is inside nnel 
!         jlev, uuint,vvint,wwint are unchanged (to avoid CPU dependency)
          if(iadptive/=0) call parallel_abort('BTRACK: adp. wrong')
          time_rm2=-99._rkind
          iexit=.true.; return
        endif !iflqs1==2

        if(iflqs1==3) then 
!         Exit during iteration in quicksearch; reduce time step and retry
          if(iadptive>=5) then
!            write(errmsg,*)'BTRACK: iadptive>=5:',iadptive,0.5_rkind*dtb,trm
!            call parallel_abort(errmsg)
            !Desperate measure
            if(trm<=0) call parallel_abort('BTRACK: trm<=0 (2d)')
            time_rm2=trm
            time_rm=time_rm-(dtb-trm)
            iexit=.true.; return
          endif !iadptive
          if(trm<=0) call parallel_abort('BTRACK: trm<=0 (2a)')
          dtb=dtb-2._rkind*trm
          dtb=dtb*real(1-1.e-3,rkind) !add safety
          if(dtb<=0) call parallel_abort('BTRACK: dtb<=0 (2b)')
          iadptive=iadptive+1
          !write(12,*)'BTRACK:',iadptive,trm,dtb,dtbk
          cycle
        endif

        iadptive=0 !reset counter

!       2nd sub-step
!       Check if vel. is too small to continue
!       This avoids 0 vel case (e.g. hits the no-slip bottom)
        uuint=uuint1; vvint=vvint1; wwint=wwint1
        vmag=sqrt(uuint*uuint+vvint*vvint)
        if(vmag<=velmin_btrack.or.iflqs1==1) exit

        xt=x0-dtb*uuint
        yt=y0-dtb*vvint
        zt=z0-dtb*wwint
        nnel=nnel0; jlev=jlev0 !reset
        call quicksearch(3,idt,l_ns,ipsgb,gcor0,frame0,dtb,x0,y0,z0,nnel,jlev, &
     &xt,yt,zt,trm,iflqs1,kbpl,arco,zrat,ztmp,vis_coe,uuint,vvint,wwint, &
     &uuint1,vvint1,wwint1)

!       Check aug. exit
        if(iflqs1==2) then 
!         Exit upon entry into quicksearch 
!         Pt (xt,yt,zt) reset, which is inside nnel; jlev reset to original
!         Forfeit the 1st half step but use the new vel. (uuint etc)
          time_rm2=-99._rkind
          iexit=.true.; return
        endif !iflqs1==2

        if(iflqs1==3) then 
!         Exit during iteration in quicksearch 
!         make sure vel. (uuint etc) pointing away from nnel
          if(trm<=0) call parallel_abort('BTRACK: trm<=0 (2b)')
          time_rm2=trm
          time_rm=time_rm-(dtb-trm)
          iexit=.true.; return
        endif

!       Check if vel. is too small to continue
!       This avoids 0 vel case (e.g. hits the no-slip bottom)
        uuint=uuint1; vvint=vvint1; wwint=wwint1
        vmag=sqrt(uuint*uuint+vvint*vvint)
        if(vmag<=velmin_btrack.or.iflqs1==1) exit


!       Update time_rm
        time_rm=time_rm-dtb
        if(time_rm<=real(1.e-6,rkind)*dt) exit

        dtb=min(dtbk,time_rm) !reset
        x0=xt
        y0=yt
        z0=zt
        nnel0=nnel; jlev0=jlev
      end do !idt
!------------------------------------------------------------------------------------------------
      endif !RK2

!     Return if for element
!     Error: Kriging for wvel as well?
      if(l_ns==3) return

      if(zrat<0._rkind.or.zrat>1._rkind) then
        write(errmsg,*)'BTRACK: zrat wrong:',jlev,zrat
        call parallel_abort(errmsg)
      endif

!     Calc max/min for ELAD
!     If inter_mom/=0, sclr() will be updated below
      if(ibtrack_test==1) then
        sclr(1)=0._rkind
        do j=1,i34(nnel)
          nd=elnode(j,nnel)
          tmp=tr_nd(1,jlev,nd)*(1._rkind-zrat)+tr_nd(1,jlev-1,nd)*zrat
          sclr(1)=sclr(1)+tmp*arco(j)
        enddo !j

        sclr(2)=-huge(1._rkind) !max
        sclr(3)=huge(1._rkind) !min
        do j=1,i34(nnel)
          nd=elnode(j,nnel)
          sclr(2)=max(sclr(2),tr_nd(1,jlev,nd),tr_nd(1,jlev-1,nd))
          sclr(3)=min(sclr(3),tr_nd(1,jlev,nd),tr_nd(1,jlev-1,nd))
        enddo !j
      else !not btrack test
        sclr(1)=-huge(1._rkind) !u max
        sclr(2)=huge(1._rkind) !u min
        sclr(3)=-huge(1._rkind) !v max
        sclr(4)=huge(1._rkind) !v min
        do j=1,i34(nnel)
          nd=elnode(j,nnel)
          sclr(1)=max(sclr(1),uu2(jlev,nd),uu2(jlev-1,nd))
          sclr(2)=min(sclr(2),uu2(jlev,nd),uu2(jlev-1,nd))
          sclr(3)=max(sclr(3),vv2(jlev,nd),vv2(jlev-1,nd))
          sclr(4)=min(sclr(4),vv2(jlev,nd),vv2(jlev-1,nd))
        enddo !j

        !Interp tracers
        sclr(5:4+ntracers)=0.d0 
        do i=1,ntracers
          do j=1,i34(nnel)
            nd=elnode(j,nnel)
            !tmp=tr_nd(i,jlev,nd)*(1._rkind-zrat)+tr_nd(i,jlev-1,nd)*zrat
            !Transport uses split, so keep vertical level at original
            tmp=tr_nd(i,j0,nd)
            sclr(4+i)=sclr(4+i)+tmp*arco(j)
          enddo !j
        enddo !i=1,ntracers
      endif !ibtrack_test

!     Kriging for vel. (excluding bnd sides)
!     For ics=2, akrmat_nd is based on eframe and need to project vel. to this frame 1st
      if(ifl_bnd/=1.and.krvel(nnel)==1) then
!       Do more inter-domain btrack if necessary to make sure the ending element is resident
        if(nnel>ne) then 
          !Check final vel. for trap
          vmag=sqrt(uuint*uuint+vvint*vvint)
          if(vmag>velmin_btrack) then 
            !Nudge the final point a little; this may create variation using different # of processors
            time_rm=real(1.e-4,rkind)*dt
            iexit=.true.; return
          endif !vmag
        else !nnel resident
!         Prepare data
          ie=ie_kr(nnel) !local index for Kriging
          if(ie==0) then
            write(errmsg,*)'Out of Kriging zone:',ielg(nnel)
            call parallel_abort(errmsg)   
          endif
          npp=itier_nd(0,ie)
          do i=1,npp
            nd=itier_nd(i,ie)
            if(idry(nd)==1) then !i.c.
              uvdata(i,1)=0._rkind
              uvdata(i,2)=0._rkind
            else !wet
!              if(ics==1) then
              vxl(1,1)=uu2(jlev,nd); vxl(1,2)=uu2(jlev-1,nd)
              vxl(2,1)=vv2(jlev,nd); vxl(2,2)=vv2(jlev-1,nd)
!              else
!                call project_hvec(uu2(jlev,nd),vv2(jlev,nd),pframe(:,:,nd),eframe(:,:,nnel),vxl(1,1),vxl(2,1))
!                call project_hvec(uu2(jlev-1,nd),vv2(jlev-1,nd),pframe(:,:,nd),eframe(:,:,nnel),vxl(1,2),vxl(2,2))
!              endif !ics
              uvdata(i,1)=vxl(1,1)*(1._rkind-zrat)+vxl(1,2)*zrat
              uvdata(i,2)=vxl(2,1)*(1._rkind-zrat)+vxl(2,2)*zrat
              !For ibtrack_test only
              uvdata(i,3)=tr_nd(1,jlev,nd)*(1._rkind-zrat)+tr_nd(1,jlev-1,nd)*zrat
            endif
          enddo !all ball nodes

          do i=1,npp+3
            al_beta(i,1:3)=0._rkind
            do j=1,npp
              al_beta(i,1:3)=al_beta(i,1:3)+akrmat_nd(i,j,ie)*uvdata(j,1:3)
            enddo !j
          enddo !i

          !Proj (xt,yt,zt) to eframe of nnel
          if(ics==1) then
            xn2=xt; yn2=yt
          else
            call project_pt('l2g',xt,yt,0._rkind,gcor0,frame0,xn1,yn1,tmp)
            call project_pt('g2l',xn1,yn1,tmp,(/xctr(nnel),yctr(nnel),zctr(nnel)/), &
     &eframe(:,:,nnel),xn2,yn2,aa1)
          endif !ics

          uuint=al_beta(npp+1,1)+al_beta(npp+2,1)*xn2+al_beta(npp+3,1)*yn2
          vvint=al_beta(npp+1,2)+al_beta(npp+2,2)*xn2+al_beta(npp+3,2)*yn2
          !For ibtrack_test=1
          if(ibtrack_test==1) sclr(1)=al_beta(npp+1,3)+al_beta(npp+2,3)*xn2+al_beta(npp+3,3)*yn2
          do i=1,npp
            nd=itier_nd(i,ie)
            if(ics==1) then
              rr=sqrt((xnd(nd)-xt)*(xnd(nd)-xt)+(ynd(nd)-yt)*(ynd(nd)-yt))
            else
              call project_pt('g2l',xnd(nd),ynd(nd),znd(nd),(/xctr(nnel),yctr(nnel),zctr(nnel)/), &
     &eframe(:,:,nnel),xn3,yn3,tmp)
              rr=sqrt((xn2-xn3)*(xn2-xn3)+(yn2-yn3)*(yn2-yn3))
            endif !ics
            covar2=covar(kr_co,rr)
            uuint=uuint+al_beta(i,1)*covar2 !dir assumed to be same as frame0 (ll)
            vvint=vvint+al_beta(i,2)*covar2
            !For ibtrack_test=1
            if(ibtrack_test==1) sclr(1)=sclr(1)+al_beta(i,3)*covar2
          enddo !i

!          !Proj vel. back to frame0
!          if(ics==2) then
!            call project_hvec(uuint,vvint,eframe(:,:,nnel),frame0,uuint1,vvint1)
!            uuint=uuint1
!            vvint=vvint1
!          endif !ics
        endif !resident element
      endif !Kriging

!     nnel wet
      end subroutine btrack

!===============================================================================
!===============================================================================

      subroutine quicksearch(idx,itr,l_ns,ipsgb,gcor0,frame0,time,x0,y0,z0,nnel, &
     &jlev,xt,yt,zt,trm,nfl,kbpl,arco,zrat,ztmp,vis_coe,uuint0,vvint0,wwint0, &
     &uuint,vvint,wwint)
!
!********************************************************************************
!										*
!     Straightline search algorithm. For ics=2, this is done in the local frame of the starting pt. 
!     
!     Inputs: 
!       idx: ID identifying where quicksearch is called (info only);
!       itr: iteration # in btrack (info only);
!       l_ns: side if l_ns=2; element if l_ns=3; (info only)
!       ipsgb: global originating node or side or element index; (info only)
!       gcor0(3): global coord. of the originating pt (for ics=2);
!       frame0(3,3): frame tensor at originating pt (2nd index is axis id) (for ics=2);
!       vis_coe: weighting factor between continuous and discontinuous vel. - not used anymore
!       time: time step from (x0,y0,z0) to (xt,yt,zt);
!       x0,y0,z0:  starting pt. (x0,y0) must be inside nnel 
!                  (for ics=2, z0 is in vertical direction and (x0,y0,z0) are in frame0); 
!       nnel,jlev: starting element and level. nnel must be inside aug. domain.
!       xt,yt,zt: projected end pt; (for ics=2, zt is in vertical direction)
!       uuint0,vvint0,wwint0: vel. used in 'time' (info only). In frame0 if ics=2.
!       In addition, su2,sv2,ww2 are also used.
! 
!     Outputs:
!      nnel, jlev: end element and level. nnel must be inside aug. domain;
!      (xt,yt,zt):  the updated end pt (if so); 
!      trm: time remaining (trm<=time). trm=0 unless the path crosses aug. bnd (nfl=3);
!      nfl: a flag. nfl=1 if a bnd or dry element is hit and vel. there is small,              
!           or death trap is reached. In this case, all outputs are valid, 
!           and the time stepping in btrack is exited successfully. 
!           If nfl=2 (hit aug. bnd upon entry) or 3 (hit aug. bnd during iteration),
!           nnel, jlev, (xt,yt,zt) and trm are updated, and nnel should be inside 
!           one of the neighbor process. If nfl=2, (xt,yt,zt) are moved to (x0,y0,z0).

!      Following outputs are valid only if nfl<=1:
!      kbpl: the local bottom level at (xt,yt);
!      arco(4): shape functions of (xt,yt);
!      zrat: vertical ratio of zt (=0 when zt=ztmp(jlev));
!      ztmp(nvrt):  z-coords. at (xt,yt);
!      uuint,vvint,wwint: interpolated vel. at end pt. In frame0 if ics=2.
!********************************************************************************
!
      use schism_glbl
      use schism_msgp, only : myrank,parallel_abort
      use hydraulic_structures, only: nhtblocks
      implicit none

      integer, intent(in) :: idx,itr,l_ns,ipsgb
      real(rkind), intent(in) :: gcor0(3),frame0(3,3),time,x0,y0,z0,vis_coe,uuint0,vvint0,wwint0
      integer, intent(inout) :: nnel,jlev
      real(rkind), intent(inout) :: xt,yt,zt
      integer, intent(out) :: nfl,kbpl
      real(rkind), intent(out) :: trm,arco(4),zrat,ztmp(nvrt),uuint,vvint,wwint

      !Function
      real(rkind) :: signa1

      !Local
      integer :: jk(4)
      real(rkind) :: wild(10,2),wild2(10,2),wild3(4,2) !,xy_l(3,2)
      real(rkind) :: vxl(4,2),vyl(4,2),vzl(4,2),vxn(4),vyn(4),vzn(4),ztmp2(nvrt,4)

      integer :: nnel00,jlev00,nel,i,j,l,nel_j,jd1,jd2,iflag,it,md1,md2,lit, &
                 &isd,n1,n2,n3,kin,nd,lev,k,inside1,inside2
      real(rkind) :: xt00,yt00,zt00,xcg,ycg,zcg,xcg2,ycg2,zcg2,pathl,ar_min1, &
                     &ar_min2,xn1,yn1,zn1,xn2,yn2,zn2,ar1,ar2,xin,yin,zin,tt1, &
                     &tt2,xt2,yt2,zt2,xtmp,ytmp,eps,dist,tmp,xctr3,yctr3,vtan, &
                     &xvel,yvel,zvel,hvel,etal,dep,hmod2,uj,vj,uj1,vj1,uu,vv,uf,vf

!     Debug
!      fdb='qs_0000'
!      lfdb=len_trim(fdb)
!      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!      open(98,file=out_dir(1:len_out_dir)//fdb,status='unknown')

!     Save for debug
      xt00=xt
      yt00=yt
      nnel00=nnel
      jlev00=jlev

!     Assumptions used in this routine:
!     (1) starting element and all crossing elements are wet
!     (2) starting pt is inside nnel; start and end pts are distinct;
!     (3) the trajectory intersects sides at 1 pt only, and when the
!         end pt (xt,yt) is near a side it's nudged into element to conclude
!         the process - this won't create CPU dependency

      if(idry_e(nnel)==1) then
        write(errmsg,*)'QUICKSEARCH: Starting element is dry:',idry_e(nnel)
        call parallel_abort(errmsg)
      endif

!     Initialize (moving) starting pt
      nfl=0
      trm=time !time remaining
      nel=nnel
      xcg=x0; ycg=y0
      pathl=sqrt((xt-xcg)*(xt-xcg)+(yt-ycg)*(yt-ycg))
      if(pathl==0._rkind.or.trm==0._rkind) then
!        write(12,*)'Last QUICKSEARCH: nodes'
!        do i=1,npa
!          do k=1,nvrt
!            write(12,*)iplg(i),k,uu2(k,i),vv2(k,i)
!          enddo !k
!        enddo !i
!        write(12,*)'Sides:'
!        do i=1,nsa
!          do k=1,nvrt
!            write(12,*)islg(i),iplg(isidenode(1:2,i)),k,su2(k,i),sv2(k,i)
!          enddo !k
!        enddo !i
        write(errmsg,*)'QUICKSEARCH: Zero path',idx,itr,l_ns,ipsgb,ielg(nel),jlev, &
     &x0,y0,xt,yt,xcg,ycg,time,uuint0,vvint0,wwint0
        call parallel_abort(errmsg)
      endif

!     Check start and end pts
      if(i34(nel)==3) then 
        call area_coord(0,nel,gcor0,frame0,xcg,ycg,arco)
        ar_min1=minval(arco(1:3)) !info for debug only
        call area_coord(0,nel,gcor0,frame0,xt,yt,arco)
        ar_min2=minval(arco(1:3))
        if(ar_min2>-small1) then
          !Fix A.C. for 0/negative
          if(ar_min2<=0) call area_coord(1,nel,gcor0,frame0,xt,yt,arco)
          nnel=nel
          trm=0._rkind
          go to 400
        endif
      else !quad
        !Reproject pt in eframe of nel for ics=2
        if(ics==1) then
          xn1=xcg; yn1=ycg
          xn2=xt; yn2=yt
        else !ll
          call project_pt('l2g',xcg,ycg,0._rkind,gcor0,frame0,xcg2,ycg2,zcg2)
          call project_pt('g2l',xcg2,ycg2,zcg2,(/xctr(nel),yctr(nel),zctr(nel)/),eframe(:,:,nel),xn1,yn1,zn1)
          call project_pt('l2g',xt,yt,0._rkind,gcor0,frame0,xcg2,ycg2,zcg2)
          call project_pt('g2l',xcg2,ycg2,zcg2,(/xctr(nel),yctr(nel),zctr(nel)/),eframe(:,:,nel),xn2,yn2,zn2)
        endif !ics
        !call quad_shape(0,1,nel,xcg,ycg,inside1,arco) !info only
        call quad_shape(0,1,nel,xn1,yn1,inside1,arco) !info only
        ar_min1=minval(arco) !info for debug only
        !call quad_shape(0,2,nel,xt,yt,inside2,arco)
        call quad_shape(0,2,nel,xn2,yn2,inside2,arco)
        ar_min2=minval(arco)
        if(inside2/=0) then
          nnel=nel
          trm=0._rkind
          go to 400
        endif
      endif !i34

!     (xt,yt) not in nel, and thus (x0,y0) and (xt,yt) are distinctive
!     Find starting edge nel_j
!     Try this twice to account for underflow (e.g. inter-btrack etc),
!     the 2nd try may create CPU dependency but it occurs very rarely
      loop6: do i=1,2
        wild=0._rkind; wild2=0._rkind !initialize for debugging output
        nel_j=0
        do j=1,i34(nel) !sides
          jd1=elnode(nxq(1,j,i34(nel)),nel)
          jd2=elnode(nxq(2,j,i34(nel)),nel)
          if(ics==1) then
            xn1=xnd(jd1); yn1=ynd(jd1)
            xn2=xnd(jd2); yn2=ynd(jd2)
          else !lat/lon
            call project_pt('g2l',xnd(jd1),ynd(jd1),znd(jd1),gcor0,frame0,xn1,yn1,zn1)
            call project_pt('g2l',xnd(jd2),ynd(jd2),znd(jd2),gcor0,frame0,xn2,yn2,zn2)
          endif !ics
          wild3(j,1)=xn1; wild3(j,2)=yn1 !save for computing centroid and nudging later
          ar1=signa1(xcg,xn1,xt,ycg,yn1,yt)    
          ar2=signa1(xcg,xt,xn2,ycg,yt,yn2)    
          wild2(j,1)=ar1; wild2(j,2)=ar2
          if(ar1>0._rkind.and.ar2>0._rkind) then
            call intersect2(xcg,xt,xn1,xn2,ycg,yt,yn1,yn2,iflag,xin,yin,tt1,tt2)
            wild(j,1)=tt1; wild(j,2)=tt2; wild(3+j,1)=xin; wild(3+j,2)=yin
            if(iflag/=1) then
              if(ics==1) then
                xcg2=xcg; ycg2=ycg; zcg2=0._rkind; xt2=xt; yt2=yt; zt2=0._rkind
              else !lat/lon
                call project_pt('l2g',xcg,ycg,0._rkind,gcor0,frame0,xcg2,ycg2,zcg2)
                call project_pt('l2g',xt,yt,0._rkind,gcor0,frame0,xt2,yt2,zt2)
              endif !ics
              write(errmsg,*)'QUICKSEARCH: Found no intersecting edges (1):',idx,itr, &
     &ielg(nel),xcg2,ycg2,zcg2,xt2,yt2,zt2,ar_min1,ar_min2,wild(1:3,1:2),wild(4:6,1:2),ar1,ar2, &
     &xcg,ycg,xt,yt,time,trm,jlev,uuint0,vvint0,wwint0,jlev00
              call parallel_abort(errmsg)
            else !success
              nel_j=j; exit loop6
            endif
          endif !ar1>=0.and.ar2>=0
        enddo !j=1,i34

        if(nel_j==0) then
          if(ics==1) then
            xcg2=xcg; ycg2=ycg; zcg2=0._rkind; xt2=xt; yt2=yt; zt2=0._rkind
          else !lat/lon
            call project_pt('l2g',xcg,ycg,0._rkind,gcor0,frame0,xcg2,ycg2,zcg2)
            call project_pt('l2g',xt,yt,0._rkind,gcor0,frame0,xt2,yt2,zt2)
          endif !ics

          if(i==1) then !1st try
            write(12,*)'QUICKSEARCH: no intersecting edge; start ID (node/side/elem)=', &
     &l_ns,'; start gb. node/side/elem #=',ipsgb,'; start level=',jlev,'; current elem=',ielg(nel), &
     &'; cg (local) coord.=',xcg2,ycg2,zcg2,'; end coord.=',xt2,yt2,zt2, &
     &'; signed areas (cg,1,t)@ nodes followed by (cg,t,2)@ nodes=',wild2(1:i34(nel),1:2), &
     &'; xcg,ycg,xt,yt=',xcg,ycg,xt,yt, &
     &'; time step from cg to t=',time,'; time remaining=',trm, &
     &'; min. area coord. for cg, t=',ar_min1,ar_min2,'; input vel=',uuint0,vvint0,wwint0, &
     &idx,itr,jlev00
            !Nudge (xcg,ycg) to off centroid (to escape some tricky underflow)
            if(i34(nel)==3) then
              wild(1,1)=real(1./3.-1.12e-2,rkind); wild(2,1)=real(1./3.-1.09e-2,rkind); wild(3,1)=1._rkind-wild(1,1)-wild(2,1) !A.C.
            else
              wild(1,1)=real(0.25-1.12e-2,rkind); wild(2,1)=real(0.25-1.09e-2,rkind); 
              wild(3,1)=real(0.25+0.937e-2,rkind); wild(4,1)=1._rkind-sum(wild(1:3,1))
            endif !i34
            xtmp=dot_product(wild(1:i34(nel),1),wild3(1:i34(nel),1))
            ytmp=dot_product(wild(1:i34(nel),1),wild3(1:i34(nel),2))
            eps=real(1.019e-2,rkind)
            xcg=(1._rkind-eps)*xcg+eps*xtmp
            ycg=(1._rkind-eps)*ycg+eps*ytmp

            if(i34(nel)==3) then
              call area_coord(0,nel,gcor0,frame0,xcg,ycg,arco)
            else
              if(ics==1) call quad_shape(0,3,nel,xcg,ycg,inside1,arco)
            endif !i34
            ar_min1=minval(arco(1:i34(nel))) !info for debug only (undefined if ics=2 and quads)
          else !i=2; out of luck
            write(errmsg,*)'QUICKSEARCH: no intersecting edge; start ID (node/side/elem)=', &
     &l_ns,'; start gb. node/side/elem #=',ipsgb,'; start level=',jlev,'; current elem=',ielg(nel), &
     &'; cg (local) coord.=',xcg2,ycg2,zcg2,'; end coord.=',xt2,yt2,zt2, &
     &'; signed areas (cg,1,t)@ nodes followed by (cg,t,2)@ nodes=',wild2(1:i34(nel),1:2), &
     &'; xcg,ycg,xt,yt=',xcg,ycg,xt,yt, &
     &'; time step from cg to t=',time,'; time remaining=',trm, &
     &'; min. area coord. for cg, t=',ar_min1,ar_min2,'; input vel=',uuint0,vvint0,wwint0, &
     &idx,itr,jlev00
            call parallel_abort(errmsg)
          endif !i
        endif !nel_j=0
      enddo loop6 !i: 2 tries

!     Check aug. exit
!     nnel, jlev, and trm are unchanged; (xt,yt,zt) moved to (x0,y0,z0)
!     to be ready for inter-subdomain tracking
      if(ic3(nel_j,nel)<0) then
        xt=x0; yt=y0; zt=z0; nnel=nel
        nfl=2; return
      endif

      zin=z0 !intialize
      it=0
      loop4: do
!----------------------------------------------------------------------------------------
      it=it+1

!     Exit loop if death trap is reached
      if(it>1000) then
!        if(ifort12(3)==0) then
!          ifort12(3)=1
        write(12,*)'QUICKSEARCH: Death trap reached'
!        endif
        nfl=1
        xt=xin
        yt=yin
        zt=zin
        nnel=nel
        trm=0._rkind !min(trm,time)
        exit loop4
      endif
      md1=elnode(nxq(1,nel_j,i34(nel)),nel)
      md2=elnode(nxq(2,nel_j,i34(nel)),nel)
      
!     Compute z position 
      dist=sqrt((xin-xt)*(xin-xt)+(yin-yt)*(yin-yt))
      tmp=min(1._rkind,dist/pathl)
      zin=zt-tmp*(zt-zin)
      trm=trm*tmp !time remaining
      pathl=dist !sqrt((xin-xt)**2+(yin-yt)**2)
      if(dist==0._rkind.or.trm==0._rkind) then
!        write(errmsg,*)'QUICKSEARCH: end pt on side:',idx,itr,l_ns,ipsgb,dist,it, &
!     &ielg(nnel00),ielg(nel),zin,time,x0,y0,xt00,yt00,xin,yin,xt,yt, &
!     &uuint0,vvint0,wwint0,ar_min2,jlev00,vis_coe
!        call parallel_abort(errmsg)

        if(i34(nel)==3) then
          call area_coord(1,nel,gcor0,frame0,xt,yt,arco) !'1' - fix A.C.
        else
          !Reproject pt in eframe of nel for ics=2
          if(ics==1) then
            xn2=xt; yn2=yt
          else !ll
            call project_pt('l2g',xt,yt,0._rkind,gcor0,frame0,xcg2,ycg2,zcg2)
            call project_pt('g2l',xcg2,ycg2,zcg2,(/xctr(nel),yctr(nel),zctr(nel)/),eframe(:,:,nel),xn2,yn2,zn2)
          endif !ics

          !call quad_shape(1,4,nel,xt,yt,inside2,arco)
          call quad_shape(1,4,nel,xn2,yn2,inside2,arco)
        endif
        nnel=nel
        trm=0._rkind
        exit loop4
      endif

!     Check for aug. exit
      if(ic3(nel_j,nel)<0) then
!       nnel is the last element inside aug. domain
!       xt,yt,zt, jlev, and trm are updated.
!       IMPORTANT: as new (xt,yt) is right on a side, make sure
!       that the vel. is pointing away from nnel!
        nfl=3
        xt=xin
        yt=yin
        zt=zin
        nnel=nel
        trm=min(trm,time) !>0
        nnel=nel
        return
      endif

!     Next element is inside aug. domain
      lit=0 !flag
!     For horizontal exit and dry elements, compute tangential vel.,
!     update target (xt,yt,zt) and continue.
!     max() added for bound check
      if(ic3(nel_j,nel)==0.or.idry_e(max(1,ic3(nel_j,nel)))==1) lit=1
      if(ihydraulics/=0.and.nhtblocks>0) then
        if(isblock_sd(1,elside(nel_j,nel))>0) lit=1 !active block
      endif
   
      if(lit==1) then
        isd=elside(nel_j,nel)
        if(isidenode(1,isd)+isidenode(2,isd)/=md1+md2) then
          write(errmsg,*)'QUICKSEARCH: Wrong side'
          call parallel_abort(errmsg)
        endif

!       Nudge intersect (xin,yin), and update starting pt
        eps=real(1.e-2,rkind) !100*small2
        if(ics==1) then
          xctr3=xctr(nel); yctr3=yctr(nel)
        else !lat/lon
!          eps=small2 !need more accuracy for south pole region
          call project_pt('g2l',xctr(nel),yctr(nel),zctr(nel),gcor0,frame0,xctr3,yctr3,tmp)
        endif !ics
        xin=(1._rkind-eps)*xin+eps*xctr3 !xctr(nel)
        yin=(1._rkind-eps)*yin+eps*yctr3 !yctr(nel)
        xcg=xin
        ycg=yin
   
        vtan=-su2(jlev,isd)*sny(isd)+sv2(jlev,isd)*snx(isd)
        !If open bnd is hit, optionally stop with nfl=1
        if(ibtrack_openbnd/=0.and.isbs(isd)>0) vtan=0._rkind
        xvel=-vtan*sny(isd)
        yvel=vtan*snx(isd)

        zvel=(ww2(jlev,md1)+ww2(jlev,md2))/2._rkind
        xt=xin-xvel*trm
        yt=yin-yvel*trm
        zt=zin-zvel*trm
        hvel=sqrt(xvel*xvel+yvel*yvel)
        if(hvel<=velmin_btrack) then
          nfl=1
          xt=xin
          yt=yin
          zt=zin
          nnel=nel
          trm=0._rkind
          exit loop4
        endif
        pathl=hvel*trm
      endif !abnormal cases

!     Search for nel's neighbor with edge nel_j, or in abnormal cases, the same element
      if(lit==0) nel=ic3(nel_j,nel) !next front element

!      do i=1,3
!        k1=elnode(i,nel)
!        k2=elnode(nx(i,1),nel)
!        if(ics==1) then
!          xn1=xnd(k1); yn1=ynd(k1)
!          xn2=xnd(k2); yn2=ynd(k2)
!        else !lat/lon
!          call project_pt('g2l',xnd(k1),ynd(k1),znd(k1),gcor0,frame0,xn1,yn1,tmp)
!          call project_pt('g2l',xnd(k2),ynd(k2),znd(k2),gcor0,frame0,xn2,yn2,tmp)
!        endif !ics
!        wild(i,1)=signa1(xn1,xn2,xt,yn1,yn2,yt)
!        !Save for debugging later
!        xy_l(i,1)=xn1; xy_l(i,2)=yn1
!      enddo !i
!      ar_min1=minval(wild(1:3,1))/area(nel)

      if(i34(nel)==3) then
        call area_coord(0,nel,gcor0,frame0,xt,yt,arco)
        ar_min1=minval(arco(1:3))
!        if(ar_min1==0) then
!          write(errmsg,*)'QUICKSEARCH impossible(2):',idx,itr,l_ns,ipsgb, &
!     &ielg(nel),ielg(nnel00),x0,y0,xt00,yt00,xt,yt,time,uuint0,vvint0,wwint0
!          call parallel_abort(errmsg)
!        endif !ar_min1==0

        if(ar_min1>-small1) then
          !arco will be fixed immediately outside loop4
          !if(ar_min1<=0) call area_coord(1,nel,gcor0,frame0,xt,yt,arco) !Fix
          nnel=nel
          trm=0._rkind
          exit loop4
        endif
      else !quad
        !Reproject pt in eframe of nel for ics=2
        if(ics==1) then
          xn2=xt; yn2=yt
        else !ll
          call project_pt('l2g',xt,yt,0._rkind,gcor0,frame0,xcg2,ycg2,zcg2)
          call project_pt('g2l',xcg2,ycg2,zcg2,(/xctr(nel),yctr(nel),zctr(nel)/),eframe(:,:,nel),xn2,yn2,zn2)
        endif !ics

        !call quad_shape(0,5,nel,xt,yt,inside2,arco)
        call quad_shape(0,5,nel,xn2,yn2,inside2,arco)
        ar_min1=minval(arco(1:4))
        if(inside2/=0) then
          !arco will be fixed immediately outside loop4
          !call quad_shape(1,?,nel,xt,yt,inside2,arco) !force the pt inside
          nnel=nel
          trm=0._rkind
          exit loop4
        endif
      endif !i34

!     Next intersecting edge
      wild=0._rkind; wild2=0._rkind !initialize for output
      nel_j=0
      do j=1,i34(nel)
        jd1=elnode(nxq(1,j,i34(nel)),nel)
        jd2=elnode(nxq(2,j,i34(nel)),nel)
!       For abnormal case, same side (border side) cannot be hit again
        if(jd1==md1.and.jd2==md2.or.jd2==md1.and.jd1==md2) cycle
        if(ics==1) then
          xn1=xnd(jd1); yn1=ynd(jd1)
          xn2=xnd(jd2); yn2=ynd(jd2)
        else !lat/lon
          call project_pt('g2l',xnd(jd1),ynd(jd1),znd(jd1),gcor0,frame0,xn1,yn1,tmp)
          call project_pt('g2l',xnd(jd2),ynd(jd2),znd(jd2),gcor0,frame0,xn2,yn2,tmp)
        endif !ics
        ar1=signa1(xcg,xn1,xt,ycg,yn1,yt)
        ar2=signa1(xcg,xt,xn2,ycg,yt,yn2)
        wild2(j,1)=ar1; wild2(j,2)=ar2
!        if(ar1>=0.and.ar2>=0) then
        if(ar1>0._rkind.and.ar2>0._rkind) then
          call intersect2(xcg,xt,xn1,xn2,ycg,yt,yn1,yn2,iflag,xin,yin,tt1,tt2)
          wild(j,1)=tt1; wild(j,2)=tt2; wild(3+j,1)=xin; wild(3+j,2)=yin
          if(iflag/=1) then
            if(ics==1) then
              xcg2=xcg; ycg2=ycg; zcg2=0._rkind; xt2=xt; yt2=yt; zt2=0._rkind
            else !lat/lon
              call project_pt('l2g',xcg,ycg,0._rkind,gcor0,frame0,xcg2,ycg2,zcg2)
              call project_pt('l2g',xt,yt,0._rkind,gcor0,frame0,xt2,yt2,zt2)
            endif !ics      
            write(errmsg,*)'QUICKSEARCH: Failed to find next edge (2):',lit,idx,itr,l_ns,ipsgb, &
     &xcg2,ycg2,zcg2,xt2,yt2,zt2,ielg(nel),iplg(md1),iplg(md2),ar_min1, &
     &wild(1:3,1:2),wild(4:6,1:2),ar1,ar2,xcg,ycg,xt,yt,time,trm,uuint0,vvint0,wwint0
            call parallel_abort(errmsg)
          else !success
            nel_j=j; !next front edge
            cycle loop4
          endif
        endif !ar1>=0.and.ar2>=0
      enddo !j

      if(nel_j==0) then
        if(ics==1) then
          xcg2=xcg; ycg2=ycg; zcg2=0._rkind; xt2=xt; yt2=yt; zt2=0._rkind
        else !lat/lon
          call project_pt('l2g',xcg,ycg,0._rkind,gcor0,frame0,xcg2,ycg2,zcg2)
          call project_pt('l2g',xt,yt,0._rkind,gcor0,frame0,xt2,yt2,zt2)
        endif !ics   
        write(errmsg,*)'QUICKSEARCH: no intersecting edge (2): ',idx,itr,l_ns,ipsgb,ielg(nel), &
     &xcg2,ycg2,zcg2,xt2,yt2,zt2,wild2(1:i34(nel),1:2),xcg,ycg,xt,yt,time,trm,uuint0,vvint0,wwint0,lit
        call parallel_abort(errmsg)
      endif !nel_j

!----------------------------------------------------------------------------------------
      end do loop4 

400   continue
!     No vertical exit from domain
      if(idry_e(nnel)==1) then
        write(errmsg,*)'QUICKSEARCH: Ending element is dry:',ielg(nnel)
        call parallel_abort(errmsg)
      endif

!     Compute area & sigma coord.
      if(i34(nnel)==3) then
        call area_coord(1,nnel,gcor0,frame0,xt,yt,arco)
      else
        !Reproject pt in eframe of nel for ics=2
        if(ics==1) then
          xn2=xt; yn2=yt
        else !ll
          call project_pt('l2g',xt,yt,0._rkind,gcor0,frame0,xcg2,ycg2,zcg2)
          call project_pt('g2l',xcg2,ycg2,zcg2,(/xctr(nnel),yctr(nnel),zctr(nnel)/),eframe(:,:,nnel),xn2,yn2,zn2)
        endif !ics

        !call quad_shape(1,6,nnel,xt,yt,inside2,arco)
        call quad_shape(1,6,nnel,xn2,yn2,inside2,arco)
      endif !i34

!     Local z-coord.
      do j=1,i34(nnel) !nodes
        nd=elnode(j,nnel)
        call zcoor(2,nd,jk(j),ztmp2(:,j))
      enddo !j
      kbpl=minval(jk(1:i34(nnel))) !min. bottom index

      ztmp(kbpl:nvrt)=0._rkind
      do j=1,i34(nnel)
        do k=kbpl,nvrt
          ztmp(k)=ztmp(k)+ztmp2(max(k,jk(j)),j)*arco(j)
        enddo !k
      enddo !j

#ifdef DEBUG
      do k=kbpl+1,nvrt
        !todo: assert
        !Warning: can be 0 for degenerate case
        if(ztmp(k)-ztmp(k-1)<0._rkind) then
          write(errmsg,*)'QUICKSEARCH: Inverted z-level in quicksearch:', &
     &ielg(nnel),k,ztmp(k),ztmp(k-1),kbpl,jk(:),ztmp(kbpl:nvrt),';', &
     &arco(1:i34(nnel)),(ztmp2(jk(j):nvrt,j),j=1,i34(nnel)), &
     &eta2(elnode(1:i34(nnel),nnel)),dp(elnode(1:i34(nnel),nnel))
          call parallel_abort(errmsg)
        endif
      enddo !k
#endif

      if(zt<=ztmp(kbpl)) then
        zt=ztmp(kbpl)
        zrat=1._rkind
        jlev=kbpl+1
      else if(zt>=ztmp(nvrt)) then
        zt=ztmp(nvrt)
        zrat=0._rkind
        jlev=nvrt
      else
        jlev=0
        do k=kbpl,nvrt-1
          if(zt>=ztmp(k).and.zt<=ztmp(k+1)) then 
            jlev=k+1
            exit
          endif
        enddo !k
        !todo: assert
        if(jlev==0) then
          write(errmsg,*)'QUICKSEARCH: Cannot find a vert. level:',zt,etal,dep,(ztmp(k),k=kbpl,nvrt)
          call parallel_abort(errmsg)
        endif
        zrat=(ztmp(jlev)-zt)/(ztmp(jlev)-ztmp(jlev-1))
      endif

      !todo: assert
      !if(zrat<0.or.zrat>1) then
      !  write(errmsg,*)'QUICKSEARCH: Sigma coord. wrong (4):',jlev,zrat
      !  call parallel_abort(errmsg)
      !endif

!      if(kbpl==kz) then !in pure S region
!        ss=(1-zrat)*sigma(jlev-kz+1)+zrat*sigma(jlev-kz)
!      else
!        ss=-99
!      endif
!      if(ss<sigma(jlev-1).or.ss>sigma(jlev)) then
!        write(11,*)'Sigma coord. wrong (5):',jlev,ss,sigma(jlev-1),sigma(jlev)
!        stop
!      endif

!     Interpolate vel.
      if(indvel==-1) then !interpolated hvel using P_1^NC
        !Pure triangles only
!       Interpolate in vertical 
        do j=1,3 !sides and nodes
          nd=elnode(j,nnel)
          isd=elside(j,nnel)
!          if(ics==1) then
          vxn(j)=su2(jlev,isd)*(1._rkind-zrat)+su2(jlev-1,isd)*zrat
          vyn(j)=sv2(jlev,isd)*(1._rkind-zrat)+sv2(jlev-1,isd)*zrat !side
!          else !lat/lon
!            call project_hvec(su2(jlev,isd),sv2(jlev,isd),sframe(:,:,isd),frame0,uj,vj)
!            call project_hvec(su2(jlev-1,isd),sv2(jlev-1,isd),sframe(:,:,isd),frame0,uj1,vj1)
!            vxn(j)=uj*(1-zrat)+uj1*zrat
!            vyn(j)=vj*(1-zrat)+vj1*zrat
!          endif !ics
          vzn(j)=ww2(jlev,nd)*(1._rkind-zrat)+ww2(jlev-1,nd)*zrat !node
        enddo !j=1,3

!       Interpolate in horizontal
        uuint=vxn(1)*(1._rkind-2._rkind*arco(1))+vxn(2)*(1._rkind-2._rkind*arco(2))+vxn(3)*(1._rkind-2._rkind*arco(3))
        vvint=vyn(1)*(1._rkind-2._rkind*arco(1))+vyn(2)*(1._rkind-2._rkind*arco(2))+vyn(3)*(1._rkind-2._rkind*arco(3))
        wwint=vzn(1)*arco(1)+vzn(2)*arco(2)+vzn(3)*arco(3)

      else !indvel>=0; interpolated hvel using P_1
!       No interpolate in time
        do j=1,i34(nnel) !nodes
          nd=elnode(j,nnel)
          do l=1,2 !levels
            lev=jlev+l-2
!            if(ics==1) then
!              uu=uu2(lev,nd); vv=vv2(lev,nd)
!            else !lat/lon
!              call project_hvec(uu2(lev,nd),vv2(lev,nd),pframe(:,:,nd),frame0,uu,vv)
!            endif !ics
            vxl(j,l)=uu2(lev,nd) !uu !+vis_coe*uf
            vyl(j,l)=vv2(lev,nd) !vv !+vis_coe*vf
            vzl(j,l)=ww2(lev,nd)
          enddo !l
        enddo !j

!       Interpolate in vertical 
        do j=1,i34(nnel)
          vxn(j)=vxl(j,2)*(1._rkind-zrat)+vxl(j,1)*zrat
          vyn(j)=vyl(j,2)*(1._rkind-zrat)+vyl(j,1)*zrat
          vzn(j)=vzl(j,2)*(1._rkind-zrat)+vzl(j,1)*zrat
        enddo !j

!       Interpolate in horizontal
        uuint=dot_product(vxn(1:i34(nnel)),arco(1:i34(nnel))) !vxn(1)*arco(1)+vxn(2)*arco(2)+vxn(3)*arco(3)
        vvint=dot_product(vyn(1:i34(nnel)),arco(1:i34(nnel))) !vyn(1)*arco(1)+vyn(2)*arco(2)+vyn(3)*arco(3)
        wwint=dot_product(vzn(1:i34(nnel)),arco(1:i34(nnel))) !vzn(1)*arco(1)+vzn(2)*arco(2)+vzn(3)*arco(3)
      endif !indvel

      end subroutine quicksearch

!===============================================================================
!===============================================================================

      subroutine intersect2(x1,x2,x3,x4,y1,y2,y3,y4,iflag,xin,yin,tt1,tt2)
!
!********************************************************************************
!										*
!     Program to detect if two segments (1,2) and (3,4) have common pts   	*
!     Assumption: the 4 pts are distinctive.					*
!     The eqs. for the 2 lines are: X=X1+(X2-X1)*tt1 and X=X3+(X4-X3)*tt2.	*
!     Output: iflag: 0: no intersection or colinear; 1: exactly 1 intersection.	*
!     If iflag=1, (xin,yin) is the intersection.				*
!     Modified to only check tt2, assuming 0<=tt1<=1 is already assured.
!										*
!********************************************************************************
!
      use schism_glbl, only: rkind 
      implicit none
      real(rkind), parameter :: small=0._rkind !small positive number or 0

      real(rkind), intent(in) :: x1,x2,x3,x4,y1,y2,y3,y4
      integer, intent(out) :: iflag
      real(rkind), intent(out) :: xin,yin,tt1,tt2

      real(rkind) :: delta,delta1,delta2

      tt1=-1000._rkind
      tt2=-1000._rkind
      xin=real(-1.e25,rkind); yin=xin
      iflag=0
      delta=(x2-x1)*(y3-y4)-(y2-y1)*(x3-x4)
      delta1=(x3-x1)*(y3-y4)-(y3-y1)*(x3-x4)
      delta2=(x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)

      if(delta/=0._rkind) then
        tt1=delta1/delta
        tt2=delta2/delta
        !if(tt1>=-small.and.tt1<=1+small.and.tt2>=-small.and.tt2<=1+small) then
        if(tt2>=-small.and.tt2<=1._rkind+small) then
          iflag=1
          xin=x3+(x4-x3)*tt2
          yin=y3+(y4-y3)*tt2
        endif
      endif

      end subroutine intersect2

!===============================================================================
!===============================================================================
! END ELCIRC BACKTRACKING SUBROUTINES
!===============================================================================
!===============================================================================

!dir$ attributes forceinline :: signa1
function signa1(x1,x2,x3,y1,y2,y3)
!-------------------------------------------------------------------------------
! Compute signed area formed by pts 1,2,3 (positive counter-clockwise)
!-------------------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg
  implicit none
  real(rkind) :: signa1
  real(rkind),intent(in) :: x1,x2,x3,y1,y2,y3

  signa1=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2._rkind
  
end function signa1

