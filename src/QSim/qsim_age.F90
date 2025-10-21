#include "fabm_driver.h"

module uzweijens_qsim_age
!
   use fabm_types
   use fabm_driver, only: type_base_driver, driver
   use schism_glbl, only: time_stamp,it_main, ne, dt ! , zone
!a
   implicit none
   private

   type, extends(type_base_model), public :: type_uzweijens_qsim_age
      ! Add variable identifiers and parameters here.
      type(type_state_variable_id) :: id_age_qsim
      type(type_global_dependency_id) :: id_day_of_year, id_seconds_of_day
      !type (type_horizontal_dependency_id) :: id_zone , id_tke_bott
      real(rk) :: age_rate
      contains
         ! Reference model procedures here.
         procedure :: initialize
         procedure :: do
   end type
   
contains

   subroutine initialize(self, configunit)
      class (type_uzweijens_qsim_age), intent(inout), target :: self
      integer,                       intent(in)            :: configunit
      character(len=256)                 :: message

      ! Register model parameters and variables here.
      call self%register_state_variable(self%id_age_qsim, "age_qsim", "d", "water_age")
      call driver%log_message('initialize - register_state_variable(self%id_age_qsim')
      call self%get_parameter(self%age_rate, "aging_rate","-","rate")
      call driver%log_message('initialize - get_parameter(self%age_rate')
      
      !call self%register_dependency(self%id_seconds_of_day, ## standard_variables%bottom_stress)
      call self%register_dependency(self%id_day_of_year, standard_variables%number_of_days_since_start_of_the_year)
      !call self%register_dependency(self%id_zone, standard_variables%bottom_stress)
      !call self%register_dependency(self%id_tke_bott, standard_variables%turbulent_kinetic_energy_at_soil_surface)
      
!models/iow/spm/spm.F90:330:   call self%register_dependency(self%id_taub, standard_variables%bottom_stress)
!models/uhh/dinoflag/dinoflag.F90:244:   call self%register_dependency(self%id_temp, standard_variables%temperature)

      !call self%register_horizontal_dependency(self%id_zone, 'zone', '-', 'Gebietsunterteilung')
      !call driver%log_message('register_horizontal_dependency(self%id_zone')
      !call self%register_dependency(self%id_zone, 'zone' , '-' , 'Zonierung des Gebiets')
      !call driver%log_message('initialize - do not register_horizontal_dependency(self%id_bottom_tke')
      !call self%register_dependency(self%id_sfpar,standard_variables%surface_downwelling_photosynthetic_radiative_flux)
      !call self%register_dependency(self%id_hs,      'snowthickness','m',    'snow thickness')

      write(message,*) 'initialized - uzweijens_qsim_age'
      call driver%log_message(message)

   end subroutine initialize

   ! Add model subroutines here.
   
   subroutine do(self, _ARGUMENTS_DO_)
      class (type_uzweijens_qsim_age), intent(in)   :: self
      _DECLARE_ARGUMENTS_DO_
      real(rk) :: age_qsim, alterung
      real(rk) :: day_of_year!, seconds_of_day
      character(len=256)                 :: message
      
      !write(message,*) it_main,' uzweijens_qsim_age do ',time_stamp,' sec. - dt=',dt
      !call driver%log_message(message)
      _GET_GLOBAL_(self%id_day_of_year,day_of_year)
      !write(message,*) ' uzweijens_qsim_age do day_of_year=',day_of_year!,' seconds_of_day=',seconds_of_day
      !call driver%log_message(message)

      _LOOP_BEGIN_
         _GET_(self%id_age_qsim,age_qsim)

         !_GET_HORIZONTAL_(self%id_zone,zone)
         !write(message,*) it_main,' uzweijens_qsim_age do zone=',zone,' | ',time_stamp,'sec.'
         !_GET_HORIZONTAL_(self%id_tke_bott, bottom_tke)
         !write(message,*) it_main,' uzweijens_qsim_age do bottom_tke=',bottom_tke,' | ',time_stamp,'sec.'
         !call driver%log_message(message)
         alterung=self%age_rate
         !! alterung=0.0_rk
         !! if((time_stamp.ge.86399.0_rk) .and. (time_stamp.le.86401.0_rk))alterung=1/dt !! hop concentration to 1 in one timestep
         !! if(abs(bottom_variable-7.0_rk) .gt. 0.1_rk) alterung=0.0_rk !! growth only in zone 7
         !if(bottom_variable- .le. 86400.0_rk) alterung=0.0_rk !!sohlvariable kann im volumen verwenden geht
         !alterung=0.0_rk
         !if(abs(ddr-7.0_rk).lt. 0.1_rk) alterung=self%age_rate
         _ADD_SOURCE_(self%id_age_qsim, alterung)
         !_ADD_SOURCE_(self%id_age_qsim, alterung*age_qsim)
         ! age_qsim=5.0_rk ! tut nix, wirft aber auch keinen fehler
         !_ADD_SOURCE_(self%id_age_qsim, self%age_rate)
         !_ADD_SOURCE_(self%id_age_qsim, self%p1*c*c*77.0_rk)
         
      _LOOP_END_
   end subroutine do

end module
