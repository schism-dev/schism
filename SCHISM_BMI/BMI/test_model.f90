module test_model

  type :: schism_type
      double precision :: model_start_time
      double precision :: model_end_time
      integer :: num_time_steps
      double precision :: current_model_time
      double precision :: time_step_size

      character(len=1000) :: SCHISM_dir

      integer :: iths
      integer :: ntime

      double precision :: ETA2
      double precision :: LatQ
      double precision :: SFCPRS
      double precision :: TMP2m
      double precision :: SOLDN
      double precision :: LWDN
      double precision :: UU10m
      double precision :: VV10m
      double precision :: Q
      double precision :: UU2
      double precision :: VV2

  end type schism_type

    contains

    subroutine run(model, dt)
        type(schism_type), intent(inout) :: model
        double precision, intent(in) :: dt

        if( dt == model%time_step_size) then
            model%ETA2 = model%ETA2 * 0.5
            model%Q = 2.0 * model%Q
            model%LatQ = 0
        else
            model%ETA2 = model%ETA2 * dt / model%time_step_size
            model%UU2 = 0.5 * dt / model%time_step_size
            model%VV2 = 0.5
        end if

        model%current_model_time = model%current_model_time + dt


    end subroutine run

end module test_model
