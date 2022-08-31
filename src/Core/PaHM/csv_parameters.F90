!----------------------------------------------------------------
!               M O D U L E   C S V _ P A R A M E T E R S
!----------------------------------------------------------------
!> @file csv_parameters.F90
!>
!> @brief
!>   Various parameters.
!>
!> @details
!>   
!>
!> @author Jacob Williams
!> @copyright License BSD
!----------------------------------------------------------------

MODULE csv_parameters

  USE PaHM_Sizes, ONLY : WP, IP

  PRIVATE

  ! maximum string length of a real number
  INTEGER(IP), PARAMETER, PUBLIC      :: max_real_str_len = 27

  ! default real number format statement (for writing real values to strings and files)
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: default_real_fmt = '(E27.17E4)'

  ! maximum string length of an integer
  INTEGER(IP), PARAMETER, PUBLIC      :: max_integer_str_len = 256

  ! default integer number format statement (for writing real values to strings and files)
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: default_int_fmt  = '(I256)'
     

END MODULE csv_parameters
