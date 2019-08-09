
! This file defines serial and parallel ASSERT macros that can be used in code or in testing. 
! 1. Why macros? 
!    a. asserts evaluate to blank lines if SELFE_ASSERTIVE is not defined. No performance hit unless they are turned on actively.
!    b. this gives us line and file data, which is very nice for failed tests.
! 2. SELFE_ASSERTIVE is always enabled for tests and will soon be enabled for Debug builds. You don't need to do anything.
! 3. There are serial (SELFE_ASSERT etc.) and parallel (MPI_ASSERT etc) flavors. The latter do some aggregating of pass/fail data in parallel.
!    You do need to pick the correct one.
! 4. There are corresponding "CHECK" macros that do not cause immediate failure. 
!    a. The test will still be indentified as a failure
!    b. Your output might be a bit more attractive by suppressing the verbose mpi abort message.



#ifdef SELFE_ASSERTIVE


!############   The ASSERTS: these are all you need from this file ############################
#define xstr(s) str(s)
#define str(s) **s


#define SELFE_CHECK(expr) call selfe_assert_impl_serial(expr,"expr",__FILE__,__LINE__,.false.)
#define SELFE_CHECK_CLOSE(val0,val1,reltol) call selfe_assert_close_impl(val0,val1,reltol,__FILE__,__LINE__,.false.,.false.)
#define SELFE_CHECK_EQUAL(val0,val1) call selfe_assert_equal_impl(val0,val1,__FILE__,__LINE__,.false.,.false.)

#define MPI_CHECK(expr) call selfe_assert_impl_mpi(expr,"expr",__FILE__,__LINE__,.false.)
#define MPI_CHECK_CLOSE(expr0,expr1,reltol) call selfe_assert_close_impl(expr0,expr1,reltol,__FILE__,__LINE__,.true.,.false.)
#define MPI_CHECK_EQUAL(expr0,expr1) call selfe_assert_equal_impl(expr0,expr1,__FILE__,__LINE__,.true.,.false.)

#define SELFE_ASSERT(expr) call selfe_assert_impl_serial(expr,"expr",__FILE__,__LINE__,.true.)
#define SELFE_ASSERT_CLOSE(val0,val1,reltol) call selfe_assert_close_impl(val0,val1,reltol,__FILE__,__LINE__,.false.,.true.)
#define SELFE_ASSERT_EQUAL(val0,val1) call selfe_assert_equal_impl(val0,val1,__FILE__,__LINE__,.false.,.true.)

#define MPI_ASSERT(expr) call selfe_assert_impl_mpi(expr,"expr",__FILE__,__LINE__)
#define MPI_ASSERT_CLOSE(expr0,expr1,reltol) call selfe_assert_close_impl(expr0,expr1,reltol,__FILE__,__LINE__,.true.,.true.)
#define MPI_ASSERT_EQUAL(expr0,expr1) call selfe_assert_equal_impl(expr0,expr1,__FILE__,__LINE__,.true.,.true.)

!##########################################################################################


! The asserts depend on these:


interface selfe_assert_close_impl
  subroutine selfe_assert_close_impl_arr(val0,val1,reltol,file,line,use_mpi,required)
  character(LEN=*), intent(in) :: file
  integer,intent(in) ::  line
  real*8, intent(in),dimension(:) ::  val0,val1
  real*8, intent(in) ::  reltol
  logical, intent(in) :: use_mpi
  logical, intent(in) :: required
  end subroutine
  subroutine selfe_assert_close_impl_dbl(val0,val1,reltol,file,line,use_mpi,required)
  character(LEN=*), intent(in) :: file
  integer,intent(in) :: line
  real*8, intent(in) ::  val0,val1
  real*8, intent(in) ::  reltol
  logical, intent(in) :: use_mpi
  logical, intent(in) :: required
  end subroutine
end interface

interface selfe_assert_equal_impl
  subroutine selfe_assert_equal_impl_arr(ival0,ival1,file,line,use_mpi,required)
  character(LEN=*), intent(in) :: file
  integer,intent(in) ::  line
  integer, intent(in),dimension(:) ::  ival0,ival1
  logical, intent(in) :: use_mpi
  logical, intent(in) :: required
  end subroutine
  subroutine selfe_assert_equal_impl_int(ival0,ival1,file,line,use_mpi,required)
  character(LEN=*), intent(in) :: file
  integer,intent(in) :: line
  integer, intent(in) ::  ival0,ival1
  logical, intent(in) :: use_mpi
  logical, intent(in) :: required
  end subroutine
end interface

#else

! If we get here, we are in optimized code and asserts are turned off by virtue of being empty macro definitions

#define SELFE_CHECK(expr) 
#define SELFE_CHECK_CLOSE(val0,val1,reltol) 
#define SELFE_CHECK_EQUAL(val0,val1) 

#define MPI_CHECK(expr) 
#define MPI_CHECK_CLOSE(expr0,expr1,reltol) 
#define MPI_CHECK_EQUAL(expr0,expr1) 

#define SELFE_ASSERT(expr) 
#define SELFE_ASSERT_CLOSE(val0,val1,reltol) 
#define SELFE_ASSERT_EQUAL(val0,val1) 

#define MPI_ASSERT(expr) 
#define MPI_ASSERT_CLOSE(expr0,expr1,reltol) 
#define MPI_ASSERT_EQUAL(expr0,expr1) 

#endif

