#include "fabm_driver.h"

module uzweijens_qsim_tracerzone
!
   use fabm_types
   use qsim_standard_variables
   use fabm_driver, only: type_base_driver, driver
   use schism_glbl, only: time_stamp,it_main, dt

   implicit none
   private

   type, extends(type_base_model), public :: type_uzweijens_qsim_tracerzone
      ! Add variable identifiers and parameters here.
      type(type_state_variable_id) :: id_tracer
      type(type_bottom_state_variable_id) :: id_bottom_variable
      type(type_global_dependency_id) :: id_day_of_year, id_seconds_of_day
      type(type_horizontal_dependency_id) :: id_zone_number

      real(rk) :: tracer_zone, tracer_time
      
      !type (type_bottom_standard_variable) :: ddr = type_bottom_standard_variable(name='ddr',  units='-')

      contains
         ! Reference model procedures here.
         procedure :: initialize
         procedure :: do
         procedure :: do_bottom
   end type type_uzweijens_qsim_tracerzone
   
contains

   subroutine initialize(self, configunit)
      class (type_uzweijens_qsim_tracerzone), intent(inout), target :: self
      integer,                       intent(in)            :: configunit
      character(len=256)                 :: message

      ! Register "planktische Variable" here.
      call self%register_state_variable(self%id_tracer, "tracer", "-", "zone marker")
      call driver%log_message('initialize - register_state_variable(self%id_tracer')
      
      ! parameters
      call self%get_parameter(self%tracer_zone, "tracered_zone","-","zone to be tracered")
      call driver%log_message('initialize - get_parameter(self%tracer_zone')
      call self%get_parameter(self%tracer_time, "tracering_time","sec.","start time of tracer event")
      call driver%log_message('initialize - get_parameter(self%tracer_time')
      
      ! time
      call self%register_dependency(self%id_day_of_year, standard_variables%number_of_days_since_start_of_the_year)
      
      ! zone in qsim_variables
      call self%register_horizontal_dependency(self%id_zone_number, qsim_variables%zone_number)

      ! test zone in bottom state_variable
      call self%register_state_variable(self%id_bottom_variable, "bottom_variable", "-", "Zone number in bottom variable")
      call driver%log_message('initialize - register_state_variable(self%id_bottom_variable')

      write(message,*) 'initialized - uzweijens_qsim_tracerzone'
      call driver%log_message(message)

   end subroutine initialize

   ! Add model subroutines here.
   
   subroutine do(self, _ARGUMENTS_DO_)
      class (type_uzweijens_qsim_tracerzone), intent(in)   :: self
      _DECLARE_ARGUMENTS_DO_
      real(rk) :: tracer, zone_hop , bottom_variable, qsim_zone
      real(rk) :: day_of_year!, seconds_of_day
      character(len=256)                 :: message
      
      !write(message,*) it_main,' type_uzweijens_qsim_tracerzone do ',time_stamp,' sec. - dt=',dt
      !call driver%log_message(message)
      _GET_GLOBAL_(self%id_day_of_year,day_of_year)
      !write(message,*) ' type_uzweijens_qsim_tracerzone do day_of_year=',day_of_year!,' seconds_of_day=',seconds_of_day
      !call driver%log_message(message)

      _LOOP_BEGIN_
         _GET_(self%id_tracer,tracer)
         _GET_HORIZONTAL_(self%id_zone_number, qsim_zone)
         _GET_HORIZONTAL_(self%id_bottom_variable, bottom_variable)
         if(abs(qsim_zone-bottom_variable).gt. 0.5 )then
            write(message,*) 'qsim_zone wrong  bottom_variable,qsim_zone=',bottom_variable,qsim_zone
            call driver%log_message(message)
         endif

         zone_hop=0.0_rk
         if((time_stamp.ge.(self%tracer_time-1)) .and. (time_stamp.le.(self%tracer_time+1)))zone_hop=1/dt !! hop concentration to 1 in one timestep
         if(abs(qsim_zone-self%tracer_zone) .gt. 0.1_rk) zone_hop=0.0_rk !! hops only in tracered_zone
         _ADD_SOURCE_(self%id_tracer, zone_hop)
         
      _LOOP_END_
   end subroutine do

subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
   class (type_uzweijens_qsim_tracerzone),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
      real(rk) :: bottom_variable

   _BOTTOM_LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_bottom_variable, bottom_variable)
         _ADD_BOTTOM_SOURCE_(self%id_bottom_variable, 0.0_rk)
   _BOTTOM_LOOP_END_
end subroutine do_bottom


end module uzweijens_qsim_tracerzone
