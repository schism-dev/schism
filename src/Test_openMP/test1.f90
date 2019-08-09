!     gfortran -ffree-line-length-none -O2 -static -static-libgfortran -finit-local-zero -fopenmp -o test1 test1.f90 
!     pgf90 -O2 -mcmodel=medium -Mbounds -mp -o test1 test1.f90
!     ifort -g -mcmodel=medium  -Bstatic -assume byterecl -qopenmp -o test1 test1.f90

!Warning: best to check if the parallel loop size >  # of threads. When they r equal, strange things happen on gcc.
!Accroding to some ref, dynamic arrays (alloc) can only be shared

      module glbl
      implicit none

      integer, save :: n
      real*8, save, allocatable :: a2d(:,:)
      end module

      program main
      use glbl
      implicit none
      real*8, parameter :: a1=1
      real*8, parameter :: b1=-1
      integer :: i,j,tid,nth,nth0,omp_get_thread_num,omp_get_num_threads,icount, & !,icount(2)
     &imax,jmax,imax2,jmax2
      real*8 :: sum1,sum2,rmax1,rmax2,rmax3
      real*8, allocatable :: aa(:,:),bb(:),cc(:),dd(:,:),ee(:),ff(:),ww(:,:)
      real*8 :: tau(3),tmp1

!     Find # of threads
!$OMP parallel default(shared)
!$OMP master
!$      nth0=omp_get_num_threads()
!$OMP end master
!$OMP end parallel

      n=5
      allocate(aa(3,n),bb(n),cc(3),dd(n,n),ee(n),ff(n*n),a2d(n,10),ww(3,n))
      bb=-1.
      do j=1,10
        do i=1,n
          a2d(i,j)=i-2*j
        enddo !i
      enddo !j

      sum1=0; sum2=0
!$OMP parallel if(n>nth0) default(shared) private(i,j,tid,nth,imax2,jmax2,rmax2,tau)
!$      tid=omp_get_thread_num()
!$      nth=omp_get_num_threads() !=1 if n<=nth0!!
!        if(nth/=nth0) then
!          print*, '# of threads mismatch,',nth,nth0
!          stop
!        endif
!$      print*, 'from thread ',tid,nth

!     The following can be init'ed outside parallel region also
!     Since they are shared, they become known to all threads after this
!$OMP master
      rmax1=-huge(1.d0)
      rmax3=rmax1
      icount=0
      imax=-1; jmax=-1
!no implied barriers for master construct
!$OMP end master
!$OMP barrier

!     Although the outer loop is parallel (and executed in random order), the inner loop is serial and
!     is executed in the intended order for correctness
!$OMP do
      do i=1,n
        ww(1,i)=-i !bottom
        do j=2,3
          ww(j,i)=ww(j-1,i)+1
        enddo !j
      enddo !i
!$OMP end do

!$OMP do reduction(+: sum1,sum2)
      !reduction vars cannot be an elem of an array e.g. rmax1(3)
      do i=1,n
        do j=1,3
          tau(j)=i+j
          print*, 'each thread',tid,i,j,tau(j)

          if(i<4.and.j<2) then
            aa(j,i)=bb(i)-n-i+j+a1
            sum1=sum1+aa(j,i)
          else
            aa(j,i)=i+j-b1
            sum2=sum2+aa(j,i)
          endif
        enddo !j
      enddo !i
!$OMP end do nowait

!$OMP workshare
      bb(1:n:2)=bb(1:n:2)-1; tmp1=maxval(abs(ww))
!$OMP end workshare

!$OMP do reduction(max: rmax1)
      do i=1,n
        imax2=-1; jmax2=-1;  rmax2=-huge(1.d0)
        do j=1,3
          if(i<5.and.j/=3) then
            rmax1=max(rmax1,aa(j,i))
!$OMP       atomic
!            icount(1)=icount(1)+1 <-atomic option cannot be applied to array elem!
            icount=icount+1 

            !Find loc in non-parallel dimension first
            if(aa(j,i)>rmax2) then
              rmax2=aa(j,i)
              imax2=i; jmax2=j
            endif
          endif !i<5.a
        enddo !j

!$OMP   critical
        if(i<5.and.rmax2>rmax3) then
          rmax3=rmax2
          imax=imax2; jmax=jmax2
        endif !aa
!$OMP   end critical
      enddo !i
!$OMP end do

!$OMP end parallel 

!     Only master thread in serial portion
       print*, 'max(ww)=',maxval(abs(ww)),tmp1
      do i=1,n
        print*, i,'ww=',ww(:,i)
        print*, i,'aa=',aa(:,i)
      enddo !i
      print*, 'sum1,2=',sum1,sum2,sum(aa)  
      print*, 'max, loc=',rmax1,rmax3,icount,imax,jmax

      print*, 'bb=',bb

!     Matrix-vector product
      call mxv(3,n,aa,bb,cc)

      print*, 'aa*bb=',cc

!     Another Matrix-vector product (no race condition)
!$OMP parallel default(shared) private(i,j,nth)
!$      nth=omp_get_num_threads()
!$      print*, 'nth(new)=',nth
!$OMP workshare
      ff=(/((i-j,i=1,n),j=1,n)/)
      dd(:,:)=reshape(ff,(/n,n/)) !transpose of the final matrix D for cache efficiency
      ee(:)=0
!$OMP end workshare

!$OMP master
      do i=1,n
        print*, i,'dd(,i)=',dd(:,i)
      enddo !i
!$OMP end master
!no implied barriers for master construct

!$OMP do schedule(guided)
      do i=1,n
        do j=1,n
          !Note the order of indices for dd; no race
          ee(i)=ee(i)+dd(j,i)*bb(j)
        enddo !j
      enddo !i
!$OMP end do
!$OMP end parallel

      print*, 'transpose(dd)*bb=',ee

!     Routine call
!$OMP parallel default(private) shared(n,ee)
      nth=omp_get_num_threads()
      tid=omp_get_thread_num()
!$OMP do
      do i=1,n
        print*, 'tst2:',tid,i
        call tst2(i,tid,ee)
      enddo !i
!$OMP end do
!$OMP end parallel

!     Nested loop
      bb=0
!$OMP parallel default(shared) private(j,i)
      do j=1,5
!$OMP   do
        do i=1,n
          bb(i)=bb(i)+i-j
        enddo !i
!$OMP   end do
      enddo !j
!$OMP end parallel

!     Compare with serial
      ee=0
      do j=1,5
        do i=1,n
          ee(i)=ee(i)+i-j
        enddo !i
      enddo !j

      print*, 'After nested loop:',bb
      print*, 'After nested loop (serial):',ee

      stop
      end

      subroutine tst2(i,tid,ee)
      use glbl
      implicit none
      integer, intent(in) :: i,tid
      real*8, intent(in) :: ee(n)

      integer :: i2,j

      i2=i*i
      write(99,*)'output from tst2 and thread ',tid,i,n,ee
      write(99,*)'i2=',i2 !still private
      do j=1,n
        write(99,*)j,'a2d:',a2d(j,:) !global a2d accesible, but better not to write to it!
      enddo !j

      end subroutine tst2

      subroutine mxv(m,n,aa,bb,cc)
      implicit none
      integer, intent(in) :: m,n
      real*8, intent(in) :: aa(m,n),bb(n)
      real*8, intent(out) :: cc(m)

      integer :: i,j

!$OMP parallel default(shared) private(i,j)
!$OMP workshare
      cc(:)=0
!$OMP end workshare

!array reduction - each thread takes a chunk of j and do array sum (1:m)
!$OMP do reduction(+: cc)
      do j=1,n
!        do i=1,m
        cc(1:m)=cc(1:m)+aa(1:m,j)*bb(j)
!        enddo !i
      enddo !j
!$OMP end parallel

      end subroutine mxv
