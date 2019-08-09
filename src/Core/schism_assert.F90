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


subroutine selfe_assert_equal_impl_int(ival0,ival1,file,line,mpi_test,required)
implicit none
character(LEN=*), intent(in) :: file
integer, intent(in) :: line
integer, intent(in) ::  ival0,ival1
logical, intent(in) :: mpi_test
character(len=1028) :: expression
logical             :: is_equal
logical :: required
is_equal = (ival0 == ival1)
write(expression,*) ival0, " == ", ival1
if (mpi_test)then
   call selfe_assert_impl_mpi(is_equal, trim(expression), file, line,required)
else
   call selfe_assert_impl_serial(is_equal,trim(expression), file, line,required)
endif
return
end subroutine

subroutine selfe_assert_equal_impl_arr(ival0,ival1,file,line,mpi_test,required)
implicit none
character(LEN=*), intent(in) :: file
integer, intent(in) :: line
integer, intent(in),dimension(:) ::  ival0,ival1
logical, intent(in) :: mpi_test
logical :: required
character(len=1028) :: expression
logical             :: is_equal
if (size(ival0) == size(ival1))then
    is_equal = maxval(abs(ival0 - ival1)) == 0
    write(expression,*) ival0, " == ", ival1, " (array comparison)"
else 
    is_equal = .false.
    write(expression,*) ival0, " == ", ival1, " (array size mismatch)"
endif
if (mpi_test)then
   call selfe_assert_impl_mpi(is_equal, trim(expression), file, line,required)
else
   call selfe_assert_impl_serial(is_equal,trim(expression), file, line,required)
endif
return
end subroutine



subroutine selfe_assert_close_impl_dbl(val0,val1,reltol,file,line,mpi_test,required)
implicit none
character(LEN=*), intent(in) :: file
integer,intent(in) :: line
real*8, intent(in) ::  val0,val1
real*8, intent(in) :: reltol
logical, intent(in) :: mpi_test
logical             :: is_close
logical :: required
character(len=1028) :: expression
is_close = (abs(val0 - val1) <= abs(reltol*val0))
write(expression,*) val0, " ~= ", val1, " reltol = ",reltol
if (mpi_test)then
   call selfe_assert_impl_mpi(is_close, trim(expression), file, line, required)
else 
   call selfe_assert_impl_serial(is_close, trim(expression), file, line, required)
end if
return
end subroutine


subroutine selfe_assert_close_impl_arr(val0,val1,reltol,file,line,mpi_test,required)
implicit none
character(LEN=*), intent(in) :: file
integer,intent(in) :: line
real*8, intent(in),dimension(:) ::  val0, val1
real*8, intent(in) ::  reltol
logical, intent(in) :: mpi_test
logical             :: is_close
logical :: required
character(len=1028) :: expression
real*8 :: scale = 0.d0
scale = maxval(abs(val0))
if (size(val0) == size(val1))then
   is_close = maxval(abs(val0 - val1)) <= abs(reltol*scale)
   write(expression,*) val0, " ~= ", val1, " reltol = ",reltol, " (array compare)"
else
   is_close = .false.
   write(expression,*) val0, " ~= ", val1, " reltol = ",reltol, " (array size mismatch)"
endif
if (mpi_test)then
   call selfe_assert_impl_mpi(is_close, trim(expression), file, line, required)
else 
   call selfe_assert_impl_serial(is_close, trim(expression), file, line, required)
end if
return
end subroutine




subroutine selfe_assert_impl_mpi(val,expression,file,line,required)
use schism_msgp
!use mpi
implicit none
include 'mpif.h'
logical val
logical is_init
logical, intent(in) :: required
character(LEN=*), intent(in) :: expression
character(LEN=*), intent(in) :: file
integer,intent(in) :: line
character(LEN=1028) :: errmsg
integer :: ierror
integer :: my_code = 0
integer :: num_fails = 0
integer :: min_failed_rank = 10000
integer :: my_rank = -1
integer :: comm_size = 0



call MPI_initialized(is_init, ierror)
if (.not. is_init .or. (ierror /= 0))then
   if (val)then
     errmsg = "Testing in parallel mode but MPI never initialize."
   else
     errmsg = "Testing in parallel mode but MPI never initialize."
   end if
   call report_fail(expression,file,line,errmsg)
   stop 1
end if

call mpi_comm_rank(comm, my_rank, ierror)
if (val) then 
  my_code=0
else
  my_code=1
end if


call MPI_allreduce(my_code,num_fails,1,MPI_INTEGER,MPI_SUM,comm,ierror)
if (num_fails == 0)then
  ! no one failed
  return
end if

if(val)then
  my_code = HUGE(my_rank)
else
  my_code = my_rank
end if


call MPI_allreduce(my_code,min_failed_rank,1,MPI_INTEGER,MPI_MIN,comm,ierror)


if (my_rank == min_failed_rank) then
  errmsg = ""
  write(errmsg,'(a,i5,8x,a,i5)') "Num processors failed: ", num_fails, "Lowest MPI rank that failed: ", min_failed_rank
  call report_fail(expression, file, line, trim(errmsg))
end if

!print*,"Required: ",required
if(required) then
  call mpi_finalize(ierror)   !(comm, -2, ierror)
  stop 2
end if



return
end subroutine


subroutine selfe_assert_impl_serial(val,expression,file,line,required)
use schism_msgp
!use mpi
implicit none
include 'mpif.h'
logical val
logical is_init
logical :: required
character(LEN=1028) :: errmsg
character(LEN=*) :: expression
character(LEN=*) :: file
integer :: line

if ( .not. val) then
    call report_fail(expression, file, line, "Serial failure")
    if (required) stop 2
end if

return
end subroutine

subroutine report_fail(expression,file,line,message)
character(LEN=*) expression, message,file
integer :: lin
logical :: required
write(*,'(a,a,a,i6,a)'), "ASSERTION FAILED: ",trim(file),"(",line,"):"
write(*,'(a)') trim(expression)
write(*,'(a)') trim(message)
write(*,*) "****"
return
end subroutine
