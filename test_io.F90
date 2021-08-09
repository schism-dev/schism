! mpif90 -cpp -O2 -mcmodel=medium -assume byterecl -g -traceback -o test_io test_io.F90
  program test_io
  implicit none
  include 'mpif.h'

  integer,parameter :: nvrt=40
  integer,parameter :: np=1000
!  integer,parameter :: mnp=10000000 !global combined max

  character(len=10) :: errchar
  integer :: ierr,comm_schism,comm2,comm,nproc_schism,myrank_schism,task_id,i,j,k,m,it 
  integer :: comm_scribe,nscribes,nproc_compute,nproc,myrank,np_gb
  integer :: nsteps,nspool,irank,noutvars,itag
  integer,allocatable :: srqst(:),sstat(:,:)
  integer,allocatable :: rrqst(:),rstat(:,:)

  real*8 :: dt,tmp,tmp_gb,tmp1
  real*8 :: work(1000), ar1(nvrt,np),ar2(nvrt,np)
  real*8, allocatable :: ar1_gb(:,:,:),ar2_gb(:,:,:)

  CALL MPI_Init(ierr)

  call mpi_comm_dup(MPI_COMM_WORLD,comm_schism,ierr)
  CALL MPI_Comm_size(comm_schism, nproc_schism, ierr)
  CALL MPI_Comm_rank(comm_schism, myrank_schism,ierr)

!  print*, 'hello from rank ',myrank_world,nproc_world

  nscribes=3
  nproc_compute=nproc_schism-nscribes
  if(myrank_schism<nproc_schism-nscribes) then !compute ranks
    task_id=1
  else !IO ranks
    task_id=2
  endif


  !Use original rank as key to order the new ranks
  CALL MPI_Comm_split(comm_schism,task_id,myrank_schism,comm2,ierr)
  CALL MPI_Comm_size(comm2, nproc, ierr)
  CALL MPI_Comm_rank(comm2, myrank,ierr)

  !comm world IDs are unreliable
  print*, 'Ranks:',myrank_schism,nproc_schism,myrank,nproc,comm2,task_id
  call mpi_barrier(comm_schism,ierr)


  !if(myrank_schism<nproc_schism-nscribes) then !compute ranks
  if(task_id==1) then !compute ranks
    comm=comm2
    if(nproc_compute/=nproc) call mpi_abort(comm_schism,'nproc mismatch',ierr)
    if(myrank+1>nproc_compute) call mpi_abort(comm_schism,'>nproc_compute',ierr)

    allocate(srqst(nscribes),sstat(MPI_STATUS_SIZE,nscribes),stat=i)
    if(i/=0) call mpi_abort(comm_schism,'alloc(1)',ierr)

    !Time step info
    dt=100.
    nsteps=4
    nspool=2
    noutvars=2
    if(noutvars>nscribes) call mpi_abort(comm_schism,'>nscribes',ierr)

    !Send to IO scribes
    if(myrank_schism==nproc_schism-nscribes-1) then
      do i=1,nscribes
        call mpi_send(dt,1,MPI_REAL8,nproc_schism-i,10,comm_schism,ierr)
        call mpi_send(nsteps,1,MPI_INTEGER,nproc_schism-i,12,comm_schism,ierr)
        call mpi_send(nspool,1,MPI_INTEGER,nproc_schism-i,13,comm_schism,ierr)
        call mpi_send(noutvars,1,MPI_INTEGER,nproc_schism-i,15,comm_schism,ierr)
      enddo !i
    endif
   
    do it=1,nsteps
      print*, 'Doing step ',it,myrank,comm,myrank_schism,nproc
      do i=1,np
        do k=1,nvrt
          ar1(k,i)=k+i+myrank+it
          ar2(k,i)=ar1(k,i)*0.5d0
        enddo !k
      enddo !i

      call mpi_barrier(comm,ierr)
      
!      tmp1=0.d0
!      do i=0,nproc-1
!        tmp1=tmp1+i**0.5d0
!      enddo !i
!      tmp=myrank**0.5d0

      if(mod(it,nspool)==0) then
        do j=1,noutvars
          print*, 'Sending to rank',nproc_schism-j,' from rank:',myrank_schism,it
          if(j==1) then
            !call mpi_send(ar1,np*nvrt,MPI_REAL8,nproc_schism-j,100+j,comm_schism,ierr)
            call mpi_isend(ar1,np*nvrt,MPI_REAL8,nproc_schism-j,100+j,comm_schism,srqst(j),ierr)
          else
            !call mpi_send(ar2,np*nvrt,MPI_REAL8,nproc_schism-j,100+j,comm_schism,ierr)
            call mpi_isend(ar2,np*nvrt,MPI_REAL8,nproc_schism-j,100+j,comm_schism,srqst(j),ierr)
          endif
          call mpi_wait(srqst(j),MPI_STATUS_IGNORE,ierr)
          !print*, 'Sent to rank',nproc_schism-j,' from rank:',myrank_schism,it
        enddo !j
        !call mpi_waitall(srqst(1:noutvars),sstat(:,1:noutvars),ierr)
        if(ierr/=MPI_SUCCESS) call mpi_abort(comm_schism,'send error',ierr)
      endif !mod()
    enddo !it
    print*, 'Done computing...',myrank
  
  else !IO ranks
    comm_scribe=comm2
    
    allocate(ar1_gb(nvrt,np,nproc_compute),ar2_gb(nvrt,np,nproc_compute))
    allocate(rrqst(nproc_compute),rstat(MPI_STATUS_SIZE,nproc_compute),stat=i)
    if(i/=0) call mpi_abort(comm_schism,'alloc(1)',ierr)

    !Get basic info
    call mpi_recv(dt,1,MPI_REAL8,nproc_schism-nscribes-1,10,comm_schism,rrqst(1),ierr)
    call mpi_recv(nsteps,1,MPI_INTEGER,nproc_schism-nscribes-1,12,comm_schism,rrqst(1),ierr)
    call mpi_recv(nspool,1,MPI_INTEGER,nproc_schism-nscribes-1,13,comm_schism,rrqst(1),ierr)
    call mpi_recv(noutvars,1,MPI_INTEGER,nproc_schism-nscribes-1,15,comm_schism,rrqst(1),ierr)

    call mpi_barrier(comm_scribe,ierr)
    print*, 'Scribe ',myrank,comm_scribe,myrank_schism
    print*, 'Scribe, basic info:',dt,nsteps,nspool,noutvars,nproc_compute

    do it=1,nsteps
!--------------------------------------------------------------------------
    if(mod(it,nspool)/=0) cycle

    print*, 'Scribe start recv...',it,myrank

    do j=1,noutvars
      if(myrank_schism==nproc_schism-j) then
        do i=1,nproc_compute
          call mpi_irecv(ar1_gb(:,:,i),np*nvrt,MPI_REAL8,i-1,100+j,comm_schism,rrqst(i),ierr)
          !call mpi_recv(ar1_gb(:,:,i),np*nvrt,MPI_REAL8,i-1,100+j,comm_schism,rrqst(i),ierr)
!        write(errchar,'(i10)')i
!        if(ierr/=MPI_SUCCESS) call mpi_abort(comm_schism,errchar,ierr)
        enddo !i
      else
        rrqst(:)=MPI_REQUEST_NULL
      endif

      call mpi_waitall(nproc_compute,rrqst,rstat,ierr)
      if(ierr/=MPI_SUCCESS) call mpi_abort(comm_schism,'receive error',ierr)

      print*, 'Scribe recv var:',it,myrank
    enddo !j
    print*, 'Scribe done recv var:',it,myrank

    if(myrank_schism==nproc_schism-1) then
      do i=1,nproc_compute
        write(98,*)'Time=',it*dt,'Compute rank ',i-1,ar1_gb(:,:,i)
        do m=1,np
          do k=1,nvrt
            if(abs(ar1_gb(k,m,i)-k-m-(i-1)-it)>1.d-5) write(97,*)'Mismatch:',i,m,k,it,ar1_gb(k,m,i),k+m+(i-1)+it
          enddo !k
        enddo !m
      enddo !i
    else if(myrank_schism==nproc_schism-2) then
      do i=1,nproc_compute
        write(99,*)'Time=',it*dt,'Compute rank ',i-1,ar1_gb(:,:,i)
        do m=1,np
          do k=1,nvrt
            if(abs(ar1_gb(k,m,i)*2.d0-k-m-(i-1)-it)>1.d-5) write(96,*)'Mismatch:',i,m,k,it,ar1_gb(k,m,i)*2
          enddo !k
        enddo !m
      enddo !i
    endif !myrank_schism

    print*, 'Scribe done writing ',it,myrank

!    do i=1,nproc_compute
!      call mpi_irecv(ar2_gb(:,:,i),np*nvrt,MPI_REAL8,i-1,14,comm_schism,rrqst(i),ierr)
!      write(errchar,'(i10)')i
!      if(ierr/=MPI_SUCCESS) call mpi_abort(comm_schism,errchar,ierr)
!    enddo !i
!    call mpi_waitall(nproc_compute,rrqst,rstat,ierr)
!    if(ierr/=MPI_SUCCESS) call mpi_abort(comm_schism,'receive error2',ierr)

!    do i=1,nproc_compute
      !Combine
!      do j=1,np
!        np_gb=j+(nproc_compute-1)*np !global node #
!        if(np_gb>mnp) call mpi_abort(comm_schism,'>mnp',ierr)
!        (:,np_gb)=ar1(:,j)

        !Debug
!        write(98,*)'Time=',it*dt,'Compute rank ',i-1,ar1_gb(:,:,i)
!        write(99,*)'Time=',it*dt,'Compute rank ',i-1,ar2_gb(:,:,i)
!      enddo !j
!    enddo !i
!--------------------------------------------------------------------------
    enddo !it
  endif



!  call MPI_Comm_free(comm_schism,ierr)
!  call MPI_Comm_free(comm2,ierr)
!  call MPI_Comm_free(comm,ierr)
!  call MPI_Comm_free(comm_scribe,ierr)

  call mpi_barrier(comm_schism,ierr)
  print*, 'Out of main time loop ',myrank_schism
  call mpi_finalize(ierr)
  if(ierr/=MPI_SUCCESS) call mpi_abort(comm_schism,'failed to finalize',ierr)
  stop

  print*, 'why am I here?'
 
  end program 
