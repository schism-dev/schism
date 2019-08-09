!!! Style comments:
Add this block, but it can be empty. License text will be automatically inserted here for distributions
!<license>
!</license>


!> Example module (in the fortran sense) containing global variables. Might also define functions to allocate them. [please expand for style comments]
!>  - Module names are lower_case_underscore with no particular designation as a module
!>  - First part of doxygen comment is prominently displayed, so choose it carefully.
!>  - One module per file preferred as it makes the build system more robust. Ignored here to keep the style guide one file.
!>  - The ingroup markup creates doxygen comments at the application module level (no relation to fortran module -- more like package)
!>  - Name collisions can occur between module and subroutine names. Avoid this by naming modules
!>     without a verb describing the action, e.g. "module gradient" and "subroutine calc_gradient"
!>  - Note that these documents use "!> " because they are above the subroutine.
!>@ingroup style_guide
module mostly_variables
    use stm_precision
    integer :: ncell  !< a module level variable documented to the right of the variable using "!<"
    
    !> Module variable documented above the variable for room using "!>"
    !> dimensions (ncell, nvar)
    real(stm_real),save,allocatable :: mass(:,:)
        
end module    


!> Example module containing routines for calculating things
!>@ingroup style_guide
module only_subroutines
contains

!> Short description of what subroutine does -- remember the first part is most visible [expand for style comments]
!>    - prefer subroutine to function
!>    - subroutine names should be lower_case_underscore and readable
!>    - variable names should be lower_case_underscore
!>    - prefer clear, slightly longer variable and routine names to "obvious" comments
!>    - output arguments should be listed before input
!>    - argument intent should be listed and args should be documented
!>    - all variables are implicit none
!>    - avoid passing anything but the most global data to subroutines by "use" statements. Passing things behind the scenes makes them
!>       untestable in isolation. 
!>    - indentat loops and conditionals two spaces
!>    - note that the !> is for comments above a subroutine and !< is for comments to the right of an argument
!>    - initialize floating point variables to LARGEREALVAL, which is a large number (1.234567d8)
!>       - do not set double precision variables to single precision constants, ie use 123.d0 not 123.0
!>       - this is more than a style issue, failure is the biggest reason debug != release
!>    - we have typesafe constants that are more fun to read: quarter, half, zero, one, ... ten and pi, gravity. 
!>    - end routines with "end subroutine" not "end".
subroutine calc_something(grad_lo,grad_hi,grad_center,vals,ncell,nvar)
use stm_precision

implicit none

!---- args
integer,intent(in)  :: ncell                          !< Number of cells
integer,intent(in)  :: nvar                           !< Number of variables
real(stm_real),intent(in)  :: vals(ncell,nvar)        !< Data to be differenced
real(stm_real),intent(out) :: grad_lo(ncell,nvar)     !< Difference on lo side, LARGEREAL in first index
real(stm_real),intent(out) :: grad_hi(ncell,nvar)     !< Difference on hi side (n+1) minus (n) LARGEREAL for last index
real(stm_real),intent(out) :: grad_center(ncell,nvar) !< Dentered diff, LARGEREAL for undefined boundary cells

!----local
integer :: ivar
S
do ivar = 1, nvar
  grad_center(2:(ncell-1),ivar) = (vals(3:ncell,ivar) - vals(1:(ncell-2),ivar))/two
  grad_center(1,ivar)=LARGEREAL
  grad_center(ncell,ivar)=LARGEREAL
  grad_hi(1:(ncell-1),ivar) = (vals(2:ncell,ivar) - vals(1:(ncell-1),ivar))
  grad_hi(ncell,ivar)=LARGEREAL
  grad_lo(2:ncell,ivar)=grad_hi(1:(ncell-1),ivar)
  grad_lo(1,ivar)=LARGEREAL
end do

return
end subroutine

end module

