

program test_hydraulic_structures
use hydraulic_structures
implicit none
#include "schism_assert.inc"


integer istruct
type(hydraulic_structure),pointer:: struct 
real(rkind),parameter :: epsilon = 5.d-16
integer, parameter    :: th_unit = 77
real(rkind)           :: tnext = 0.0d0
real(rkind)           :: tm = -1.0d0

! Creat the test harness
call load_structures('hydraulic_structure_test.in')
SELFE_ASSERT_EQUAL(nhtblocks,5)
SELFE_ASSERT_CLOSE(block_nudge,0.01d0,epsilon)

! Some structures of each of the five types
! HYDTRANSFER = 1  
! ORIFICE     = 2  
! WEIR        = 3  
! PIPE        = 4  
! RADIAL      = 5
do istruct = 1,nhtblocks
  struct => structures(istruct)
  SELFE_ASSERT_EQUAL(struct%struct_type, istruct)
  SELFE_ASSERT_EQUAL(struct%nduplicate,mod(istruct,2)+1)
  if (istruct == PIPE) then
    SELFE_ASSERT_CLOSE(struct%width,2.0d0,epsilon)
    SELFE_ASSERT(struct%height > 1.d6)
  else if(istruct == WEIR) then
    SELFE_ASSERT_CLOSE(struct%width,10.1d0,epsilon)
    SELFE_ASSERT(struct%height > 1.d16)
  else if (istruct /= HYDTRANSFER) then
    SELFE_ASSERT_CLOSE(struct%width,10.1d0,epsilon)
    SELFE_ASSERT_CLOSE(struct%height, 4.2d0,epsilon)    
  end if

  if (istruct == HYDTRANSFER)then
    SELFE_ASSERT_CLOSE(struct%prescribed_flow,5.d0,epsilon)
  else
    SELFE_ASSERT(struct%prescribed_flow==0.d0)
    SELFE_ASSERT_CLOSE(struct%elev,-1.d0,epsilon)
    SELFE_ASSERT_CLOSE(struct%op_down,0.5d0,epsilon)
    SELFE_ASSERT_CLOSE(struct%op_up,1.2d0,epsilon)
    SELFE_ASSERT_CLOSE(struct%coef,0.75d0,epsilon)
  endif
end do

open(th_unit,file='struct_one.th',status='old')
call irreg_time_history_advance(th_unit, 0.d0, tnext)
read(th_unit,*)tm
SELFE_ASSERT_CLOSE(tm,0.d0,epsilon)
SELFE_ASSERT_CLOSE(tnext,150.d0,epsilon)
call irreg_time_history_advance(th_unit, 120.d0, tnext)
SELFE_ASSERT_CLOSE(tnext,150.d0,epsilon)
call irreg_time_history_advance(th_unit, 150.d0, tnext)
SELFE_ASSERT_CLOSE(tnext,240.d0,epsilon)
call irreg_time_history_advance(th_unit, 240.d0, tnext)
SELFE_ASSERT_CLOSE(tnext,241.d0,epsilon)
call irreg_time_history_advance(th_unit, 250.d0, tnext)
SELFE_ASSERT_CLOSE(tnext,360.d0,epsilon)
call irreg_time_history_advance(th_unit, 480.d0, tnext)
SELFE_ASSERT_CLOSE(tnext,481.d0,epsilon)
close(th_unit)

open(th_unit,file='struct_one.th',status='old')
call irreg_time_history_advance(th_unit, 481.d0, tnext)
SELFE_ASSERT(tnext > 1.d10)
read(th_unit,*)tm
SELFE_ASSERT_CLOSE(tm, 481.d0,epsilon)
close(th_unit)

open(th_unit,file='struct_one.th',status='old')
call irreg_time_history_advance(th_unit, 485.d0, tnext)
SELFE_ASSERT(tnext > 1.d10)
read(th_unit,*)tm
SELFE_ASSERT_CLOSE(tm, 481.d0,epsilon)
close(th_unit)
end program



      
