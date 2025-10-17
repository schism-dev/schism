module uzweijens_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   use uzweijens_qsim_age
   ! Add use statements for new models here

   implicit none

   private

   type, extends(type_base_model_factory) :: type_factory
   contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: uzweijens_model_factory

contains

   subroutine create(self, name, model)

      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
      case ('qsim_age'); allocate (type_uzweijens_qsim_age::model)
      ! Add case statements for new models here
      end select

   end subroutine create

end module