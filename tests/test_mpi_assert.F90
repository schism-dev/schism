

!< Tests all variations where the asserts pass
program test_mpi_assert
use schism_msgp
implicit none
#include "schism_assert.inc"
real*8  :: val0, val1
integer :: int0, int1
real*8  :: arr0(4), arr1(4)
integer :: iarr0(4),iarr1(4)
real*8, parameter :: eps = 2.d-15

call parallel_init

val0 = 1.d0 + 1.d-15
val1 = 1.d0
MPI_ASSERT_CLOSE(val0,val1, eps )
MPI_CHECK_CLOSE(val0,val1, eps )


int0 = 1
int1 = int0
MPI_ASSERT_EQUAL(int0,int1)
MPI_CHECK_EQUAL(int0,int1)

arr0 = 3.7d0
arr1 = 3.7d0
arr1(3) = 3.7d0+1.d-15
MPI_ASSERT_CLOSE(arr0,arr1,eps)
MPI_CHECK_CLOSE(arr0,arr1,eps)


iarr0 = myrank
iarr1 = myrank
MPI_ASSERT_EQUAL(iarr0,iarr1)
MPI_CHECK_EQUAL(iarr0,iarr1)


MPI_ASSERT(myrank < 100 )
!MPI_CHECK( myrank < 1000 )
call parallel_finalize
end program
