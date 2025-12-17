#include "fabm_driver.h"

module uzweijens_qsim_decay
!
   use fabm_types
   use qsim_standard_variables
   use fabm_driver, only: type_base_driver, driver

   implicit none
   private

   type, extends(type_base_model), public :: type_uzweijens_qsim_decay
      ! Add variable identifiers and parameters here.
      type(type_state_variable_id) :: id_tracer
      real(rk) :: decay_rate
      contains
         ! Reference model procedures here.
         procedure :: initialize
         procedure :: do
   end type type_uzweijens_qsim_decay
   
contains

   subroutine initialize(self, configunit)
      class (type_uzweijens_qsim_decay), intent(inout), target :: self
      integer,                       intent(in)            :: configunit
      character(len=256)                 :: message

      ! Register "planktische Variable" here.
      call self%register_state_variable(self%id_tracer, "tracer", "-", "zone marker")
      call driver%log_message('initialize_decay - register_state_variable(self%id_tracer')
      ! call self%register_state_dependency(self%id_tracer, "tracer", "-", "zone marker")
      ! call driver%log_message('initialize_decay - register_state_dependency(self%id_tracer')

      ! parameters
      call self%get_parameter(self%decay_rate, "decay_rate","-","rate of tracer decay")
      call driver%log_message('initialize_decay - get_parameter(self%decay_rate')
      
      write(message,*) 'initialized - uzweijens_qsim_decay'
      call driver%log_message(message)

   end subroutine initialize

   ! Add model subroutines here.
   
   subroutine do(self, _ARGUMENTS_DO_)
      class (type_uzweijens_qsim_decay), intent(in)   :: self
      _DECLARE_ARGUMENTS_DO_
      real(rk) :: tracer
      character(len=256)                 :: message
      
      _LOOP_BEGIN_
         _GET_(self%id_tracer,tracer)
         _ADD_SOURCE_(self%id_tracer, -1.0_rk * self%decay_rate * tracer)
      _LOOP_END_
   end subroutine do

end module uzweijens_qsim_decay
