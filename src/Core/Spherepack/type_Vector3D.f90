!
!  Purpose:
!
!  Defines a class for 3-dimensional cartesian vector calculations
!
module type_Vector3D

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: assignment(=)
    public :: operator(*)
    
    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp

    type, public :: Vector3D
        ! Type components
        real(wp), public :: x = ZERO
        real(wp), public :: y = ZERO
        real(wp), public :: z = ZERO
    contains
        ! Type-bound procedures
        procedure,          private  :: add_vectors
        procedure,          private  :: subtract_vectors
        procedure,          private  :: divide_vector_by_real
        procedure,          private  :: divide_vector_by_integer
        procedure,          private  :: get_dot_product
        procedure,          private  :: assign_vector_from_integer_array
        procedure,          private  :: assign_vector_from_real_array
        procedure, nopass,  private  :: assign_real_array_from_vector
        procedure,          private  :: copy_vector
        procedure,          private  :: multiply_vector_times_real
        procedure, nopass,  private  :: multiply_real_times_vector
        procedure,          private  :: multiply_vector_times_integer
        procedure, nopass,  private  :: multiply_integer_times_vector
        procedure,          private  :: get_cross_product
        procedure,          public   :: get_norm
        ! Generic type-bound procedures
        generic, public :: operator(.dot.) => get_dot_product
        generic, public :: operator(.cross.) => get_cross_product
        generic, public :: operator(+) => add_vectors
        generic, public :: operator(-) => subtract_vectors
        generic, public :: operator(/) => &
            divide_vector_by_real, &
            divide_vector_by_integer
    end type Vector3D

    type, public :: Vector3DPointer
        ! Type components
        type(Vector3D), pointer :: ptr => null()
    end type Vector3DPointer

    ! Interface for assignment operator
    interface assignment(=)
        module procedure assign_vector_from_integer_array
        module procedure assign_vector_from_real_array
        module procedure assign_real_array_from_vector
        module procedure copy_vector
    end interface

    ! Interface for multiplication operator
    interface operator(*)
        module procedure multiply_vector_times_real
        module procedure multiply_real_times_vector
        module procedure multiply_vector_times_integer
        module procedure multiply_integer_times_vector
        module procedure get_cross_product
    end interface

contains

    subroutine assign_vector_from_integer_array(self, array)

        ! Dummy arguments
        class(Vector3D),  intent(out) :: self
        integer(ip),      intent(in)  :: array(:)

        select type(self)
            type is (Vector3D)
            self%x = real(array(1), kind=wp)
            self%y = real(array(2), kind=wp)
            self%z = real(array(3), kind=wp)
        end select

    end subroutine assign_vector_from_integer_array

    subroutine assign_vector_from_real_array(self, array)

        ! Dummy arguments
        class(Vector3D), intent(out) :: self
        real(wp),        intent(in)  :: array(:)

        select type(self)
            type is (Vector3D)
            self%x = array(1)
            self%y = array(2)
            self%z = array(3)
        end select
    end subroutine assign_vector_from_real_array

    pure subroutine assign_real_array_from_vector(array, self)

        ! Dummy arguments
        real(wp),        intent(out) :: array(:)
        class(Vector3D), intent(in)  :: self

        select type(self)
            type is (Vector3D)
            array(1) = self%x
            array(2) = self%y
            array(3) = self%z
        end select

    end subroutine assign_real_array_from_vector

    subroutine copy_vector(self, other)

        ! Dummy arguments
        class(Vector3D), intent(out) :: self
        class(Vector3D), intent(in)  :: other

        self%x = other%x
        self%y = other%y
        self%z = other%z

    end subroutine copy_vector

    pure function add_vectors(self, other) &
        result (return_value)

        ! Dummy arguments
        class(Vector3D), intent(in) :: self
        class(Vector3D), intent(in) :: other
        type(Vector3D)              :: return_value

        return_value%x = self%x + other%x
        return_value%y = self%y + other%y
        return_value%z = self%z + other%z

    end function add_vectors

    pure function subtract_vectors(self, other) &
        result (return_value)

        ! Dummy arguments
        class(Vector3D), intent(in) :: self
        class(Vector3D), intent(in) :: other
        type(Vector3D)              :: return_value

        return_value%x = self%x - other%x
        return_value%y = self%y - other%y
        return_value%z = self%z - other%z

    end function subtract_vectors

    pure function multiply_vector_times_real(self, real_var) &
        result (return_value)

        ! Dummy arguments
        class(Vector3D), intent(in) :: self
        real(wp),        intent(in) :: real_var
        type(Vector3D)              :: return_value

        return_value%x = self%x * real_var
        return_value%y = self%y * real_var
        return_value%z = self%z * real_var

    end function multiply_vector_times_real

    pure function multiply_real_times_vector(real_var, other) &
        result (return_value)

        ! Dummy arguments
        real(wp),        intent(in) :: real_var
        class(Vector3D), intent(in) :: other
        type(Vector3D)              :: return_value

        return_value%x = real_var * other%x
        return_value%y = real_var * other%y
        return_value%z = real_var * other%z

    end function multiply_real_times_vector

    pure function multiply_vector_times_integer(self, int_var) &
        result (return_value)

        ! Dummy arguments
        class(Vector3D), intent(in) :: self
        integer(ip),     intent(in) :: int_var
        type(Vector3D)              :: return_value

        ! Local variables
        real(wp) :: real_var

        real_var = real(int_var, kind=wp)

        return_value%x = self%x * real_var
        return_value%y = self%y * real_var
        return_value%z = self%z * real_var

    end function multiply_vector_times_integer

    pure function multiply_integer_times_vector(int_var, other) &
        result (return_value)

        ! Dummy arguments
        integer(ip),     intent(in) :: int_var
        class(Vector3D), intent(in) :: other
        type(Vector3D)              :: return_value

        ! Local variables
        real(wp) :: real_var

        real_var = real(int_var, kind=wp)

        return_value%x = real_var * other%x
        return_value%y = real_var * other%y
        return_value%z = real_var * other%z

    end function multiply_integer_times_vector

    pure function divide_vector_by_real(self, real_var) &
        result (return_value)

        ! Dummy arguments
        class(Vector3D), intent(in) :: self
        real(wp),        intent(in) :: real_var
        type(Vector3D)              :: return_value

        return_value%x = self%x / real_var
        return_value%y = self%y / real_var
        return_value%z = self%z / real_var

    end function divide_vector_by_real

    pure function divide_vector_by_integer(self, int_var) &
        result (return_value)

        ! Dummy arguments
        class(Vector3D), intent(in) :: self
        integer(ip),     intent(in) :: int_var
        type(Vector3D)              :: return_value

        return_value%x = self%x / int_var
        return_value%y = self%y / int_var
        return_value%z = self%z / int_var

    end function divide_vector_by_integer

    pure function get_dot_product(self, other) &
        result (return_value)

        ! Dummy arguments
        real(wp)                    :: return_value
        class(Vector3D), intent(in) :: self
        class(Vector3D), intent(in) :: other

        associate (&
            a => [self%x, self%y, self%z], &
            b => [other%x, other%y, other%z] &
           )
            return_value = dot_product(a, b)
        end associate

    end function get_dot_product

    pure function get_cross_product(self, other) &
        result (return_value)

        ! Dummy arguments
        class(Vector3D), intent(in) :: self
        class(Vector3D), intent(in) :: other
        type(Vector3D)              :: return_value

        return_value%x = self%y * other%z - self%z * other%y
        return_value%y = self%z * other%x - self%x * other%z
        return_value%z = self%x * other%y - self%y * other%x

    end function get_cross_product

    pure function get_norm(self) &
        result (return_value)

        ! Dummy arguments
        class(Vector3D), intent(in) :: self
        real(wp)                    :: return_value

        associate (array => [self%x, self%y, self%z])
            return_value = norm2(array)
        end associate

    end function get_norm

end module type_Vector3D
