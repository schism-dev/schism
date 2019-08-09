program combine_icm_output
!combine ICM model station output and dump results into one file
!
!on viz3
!ifort -O2 -Bstatic -assume byterecl -o combine_icm_output combine_icm_output.f90

  implicit none
  integer, parameter :: rkind=8
  integer :: i,j,k,m,id,it,n,n1,n2,istat
  integer :: nsgb,ne,np,nvrt,nproc,nstation,ntr,nday,iform
  integer :: ftell
  real(rkind) :: t0,time,dt
  character(len=4) :: fn
  logical :: lexist

  integer, allocatable :: ns(:),iexist(:)
  real(rkind), allocatable :: S(:,:)
  character(len=12),allocatable :: trname(:)
 
  !read number for rank 
  open(10,file='./outputs/local_to_global_0000',status='old');
  read(10,*)nsgb,ne,np,nvrt,nproc
  close(10)
  
  !read station info.
  open(11,file='cstation.in',status='old')
  read(11,*); read(11,*)nstation
  close(11)

  !read in # of tracer and nday
  write(*,*)"input: number of tracers"
  read(*,*)ntr
  write(*,*)"input: number of days (end point)"
  read(*,*)nday
  write(*,*)"input: format of result (1:ASCII; 2: binary) "
  read(*,*)iform

  !allocation
  allocate(S(ntr,nstation),ns(nproc),iexist(nproc),trname(ntr), stat=istat)
  if(istat/=0) stop 'failed in alloc. S'

  !read raw output information
  t0=1e16
  iexist=0 
  do i=1,nproc
    write(fn,'(i4.4)')i-1
    inquire(file='./outputs/cstation_'//fn//'.out',exist=lexist)
    if(lexist) then
      iexist(i)=1
      open(20+i,file='./outputs/cstation_'//fn//'.out',form='unformatted',status='old')
      read(20+i)ns(i),dt
      n1=ftell(20+i)
      read(20+i)time,id,(S(m,id),m=1,ntr) 
      n2=ftell(20+i)

      t0=min(t0,time)
      call fseek(20+i,n1-n2,1)
    endif
  enddo

  !oopen file to write
  if(iform==1) then !ascii
    open(12,file='outputs/cstation.out',status='replace')
  elseif(iform==2) then !binary
    open(12,file='outputs/cstation.out',access='stream',form='unformatted',status='replace')
  else
    stop 'unknown iform'
  endif

  !read result and dump into one cstation.out
  do while(.true.)
    if(t0>nday*86400) exit

    S=-9999
    do i=1,nproc
      if(iexist(i)==0) cycle
      do j=1,ns(i)
        n1=ftell(20+i)
        read(20+i,iostat=istat)time,id,(S(m,id),m=1,ntr) 
        if(istat/=0) then
          close(20+i)
          iexist(i)=0
        endif
        n2=ftell(20+i)
        if(time>t0) then
          call fseek(20+i,n1-n2,1)
          exit 
        endif
      enddo !j
    enddo !i

    !dump results
    if(iform==1) then
      write(12,'(f12.4,x,100000(f14.5,x))')t0/86400,((S(m,i),i=1,nstation),m=1,ntr)
    elseif(iform==2) then
      write(12)t0/86400,((S(m,i),i=1,nstation),m=1,ntr)
    endif

    
    t0=t0+dt
  enddo !while

  stop
end program combine_icm_output
