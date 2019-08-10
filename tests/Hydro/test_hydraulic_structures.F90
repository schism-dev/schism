

program test_hydraulic_structures
use hydraulic_structures
implicit none
#include "schism_assert.inc"

integer, parameter :: nstruct = 6
integer istruct
real(rkind) :: flow,flow2,downelev,upelev,frac
real(rkind) :: rhcoef,subratio
real(rkind) :: outside_calc
real(rkind) :: epsilon = 1.d-13
type(hydraulic_structure),pointer:: struct 
! use of intrinsics in initialization sqrt(2.d0*gravity) not allowed by pgi
real(rkind),parameter :: sqrt2grav = 4.428690551d0 
real(rkind), parameter   :: PART_SUBMERGE = 2.d0/3.d0
real(rkind), parameter   :: FULL_SUBMERGE = 0.8d0

! Creat the test harness
! Some structures of each of the five types
! HYDTRANSFER = 1  
! ORIFICE     = 2  
! WEIR        = 3  
! PIPE        = 4  
! RADIAL      = 5
! RADIAL_RH   = 6
allocate(structures(nstruct))
do istruct = 1,nstruct
  struct => structures(istruct)
  struct%struct_type = istruct
  struct%elev = -1.d0
  if (istruct == PIPE) then
    struct%width = 2.0d0
    struct%height = 2.d0*struct%width
  end if
  struct%width = 10.d0
  struct%height = 4.0d0
  struct%op_down = 0.5d0
  struct%op_up = 1.0d0
  struct%coef = 0.75d0
  struct%coef_lin = 0.25d0
  struct%prescribed_flow = 5.d0  
end do


upelev = 1.0d0
downelev = 0.75d0
call calc_struc_flow(1,upelev,downelev,flow)
SELFE_ASSERT_CLOSE(flow,5.0d0,1.d-15)

call calc_struc_flow(2,upelev,downelev,flow)
outside_calc =  10.d0*(upelev+1.d0)*sqrt2grav*sqrt(upelev-downelev)*0.5d0*0.75d0
SELFE_ASSERT_CLOSE(flow,outside_calc,epsilon)


!reverse is negative and double magnitude because of up_down vs op_up
call calc_struc_flow(2,downelev,upelev,flow)
SELFE_ASSERT_CLOSE(flow,-2.d0*outside_calc,epsilon)

! submerged flow where up and down elevs are above top of orifice
structures(2)%height = 1.0d0
outside_calc =  10.d0*(1.d0)*sqrt2grav*sqrt(upelev-downelev)*0.5d0*0.75d0
call calc_struc_flow(2,upelev,downelev,flow)
SELFE_ASSERT_CLOSE(flow,outside_calc,epsilon)


! still submerged upstream, not down
downelev = -1.25d0
call calc_struc_flow(2,upelev,downelev,flow)
outside_calc =  10.d0*(1.d0)*sqrt2grav*sqrt(upelev+1.d0)*0.5d0*0.75d0
SELFE_ASSERT_CLOSE(flow,outside_calc,epsilon)

! everything below the orifice invert elevation
upelev = -1.25d0
call calc_struc_flow(2,upelev,downelev,flow)
SELFE_ASSERT_CLOSE(flow,0.d0,epsilon)

! Reverse submerged flow, up elev below invert orifice flow, 2 duplicates 
upelev = 1.0d0
structures(2)%nduplicate=2
call calc_struc_flow(2,downelev,upelev,flow)
outside_calc =  -2.d0*10.d0*(1.d0)*sqrt2grav*sqrt(upelev+1.d0)*1.0d0*0.75d0
SELFE_ASSERT_CLOSE(flow,outside_calc,epsilon)


! Unsubmerged weir flow
upelev = 1.0d0
downelev = 0.75d0
structures(3)%height = huge(1.0d0)
call calc_struc_flow(3,upelev,downelev,flow)
outside_calc =  10.d0*2.d0*sqrt2grav*sqrt(upelev-downelev)*0.5d0*0.75d0
!!SELFE_ASSERT_CLOSE(flow,outside_calc,epsilon)

! Reverse unsubmerged weir flow, with op_up = 2*op_down
call calc_struc_flow(3,downelev,upelev,flow)
outside_calc =  -2.d0*outside_calc
!!!!!SELFE_ASSERT_CLOSE(flow,outside_calc,epsilon)



!pathological case -- this is just here as a reminder
structures(3)%height = 1.0d0
call calc_struc_flow(3,upelev,downelev,flow)
print*,"Upper submerged height weir flow",flow

downelev = -1.25d0
structures(3)%height = huge(1.0d0)
call calc_struc_flow(3,upelev,downelev,flow)
print*,"Low elev below invert submerged weir flow",flow
upelev = -1.25d0
call calc_struc_flow(3,upelev,downelev,flow)
print*,"Both elev below invert weir flow",flow
upelev = 1.0d0
call calc_struc_flow(3,downelev,upelev,flow)
structures(3)%nduplicate=2
print*,"Reverse submerged flow, up elev below invert weir flow, 2 duplicates",flow




print*,"********************"
upelev = 1.0d0
downelev = 0.25d0
structures(5)%elev = -1.d0
structures(5)%height = 1.d0


call calc_struc_flow(5,upelev,downelev,flow)
subratio = 1.25d0/2.d0   ! 0.625, just below threshold of 0.67
! width*flow_aperture*op*coef
outside_calc = 10.d0*1.d0*0.5d0*0.75d0*sqrt2grav*sqrt(upelev-(-1.d0))
SELFE_ASSERT_CLOSE(flow,outside_calc,epsilon)
print*,"Almost partially submerged free radial gate flow",flow

structures(5)%height = 2.d0
call calc_struc_flow(5,upelev,downelev,flow)
outside_calc = 10.d0*2.d0*0.5d0*0.75d0*sqrt2grav*sqrt(upelev-(-1.d0))
SELFE_ASSERT_CLOSE(flow,outside_calc,epsilon)
print*,"Totally Free flow radial gate flow",flow


! Test transition between partially submerged and free. Assumes transitions constants of 
! for PART_SUBMERGED and FULL_SUBMERGED, so if these depart from current hardwired values
! this test would fail
upelev = 1.0d0
downelev = 0.5999999999d0
structures(5)%height = 1.d0
call calc_struc_flow(5,upelev,downelev,flow)
print*,"Partially submerged radial gate flow",flow


upelev = 1.0d0
downelev = 0.60000000001d0
! this is right at the borderline between free and partially submerged
call calc_struc_flow(5,upelev,downelev,flow2)
print*,"Free radial gate flow",flow
SELFE_ASSERT_CLOSE(flow, flow2, 5.d-3)


upelev = 1.0d0
downelev = 0.5d0 
! submergence is 0.75d0
frac = (0.75-PART_SUBMERGE)/(FULL_SUBMERGE-PART_SUBMERGE)
call calc_struc_flow(5,upelev,downelev,flow)
! width*aperture_height*op"*coef
outside_calc = frac*10.d0*1.d0*0.5*0.75*sqrt2grav*sqrt(upelev-downelev) &
       +(1.d0-frac)*10.d0*1.d0*0.5*0.75*sqrt2grav*sqrt(upelev-downelev)*sqrt(3.d0)
SELFE_ASSERT_CLOSE(flow,outside_calc,epsilon)

! TEst transition btween partially submerged and fully
print*,"Partially submerged radial gate flow",flow
upelev = 1.0d0
downelev = 0.3333333d0
call calc_struc_flow(5,upelev,downelev,flow)
print*,"Partially submerged radial gate flow",flow

upelev = 1.0d0
downelev = 0.333333334d0
call calc_struc_flow(5,upelev,downelev,flow2)
SELFE_ASSERT_CLOSE(flow, flow2, 5.d-3)


downelev = -1.25d0
call calc_struc_flow(5,upelev,downelev,flow)
print*,"Low elev below invert unsubmerged radial flow",flow
upelev = -1.25d0
call calc_struc_flow(5,upelev,downelev,flow)
print*,"Both elev below invert radial flow",flow
upelev = 1.0d0

print*,"********************"
structures(4)%width  = 2.d0
structures(4)%height = 4.d0

upelev = 1.0d0
downelev = 0.25d0
call calc_struc_flow(4,upelev,downelev,flow)
print*,"Pipe: half submerged headwater, tail",flow

upelev = 3.0d0
call calc_struc_flow(4,upelev,downelev,flow)
print*,"Pipe: Submerged headwater, partly submerged tail flow",flow

upelev = 3.0d0
downelev = -1.25d0
call calc_struc_flow(4,upelev,downelev,flow)
print*,"Pipe: submerged head, unsubmerged tail ",flow

downelev = 3.0d0
upelev = -1.25d0
call calc_struc_flow(4,upelev,downelev,flow)
print*,"Pipe: Unsubmerged head, submerged tail ",flow

upelev = -1.25d0
downelev = -1.25d0
call calc_struc_flow(4,upelev,downelev,flow)
print*,"Pipe: both elev below invert flow",flow

SELFE_ASSERT( upelev==downelev )
!SELFE_ASSERT_CLOSE(upelev,downelev+1.d-5,1.d-15)


upelev = 5.d0
downelev = 4.5d0
call calc_struc_flow(6,upelev,downelev,flow)
! gate height is 4, relative gate height is 4/6 where 6 is 5-(-1)
rhcoef = 0.75d0 + 0.25d0*(2d0/3.d0)
! width*height_of_flow*op*coef
outside_calc =  10.d0*4.d0*5.d-1*rhcoef*sqrt2grav*sqrt(upelev-downelev)
SELFE_ASSERT_CLOSE(flow,outside_calc,epsilon)

upelev = 2.0d0
downelev = 1.5d0
call calc_struc_flow(6,upelev,downelev,flow)
! gate height is 4, relative gate height is 1.0 because of the ceiling
rhcoef = 1.d0
! width*height_of_flow*op*coef
outside_calc =  10.d0*3.d0*5.d-1*rhcoef*sqrt2grav*sqrt(upelev-downelev)
SELFE_ASSERT_CLOSE(flow,outside_calc,epsilon)



end program
